import sys,json,os,numpy as np,pandas as pd
from rdkit import Chem,RDLogger; from rdkit.Chem import rdFingerprintGenerator as G
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
RAW=sys.argv[1]; OUT=sys.argv[2]; os.makedirs(OUT,exist_ok=True)
T=json.load(open('data/targets_resolved.json')); T['ADRA1A']=[['CHEMBL229','a']]
gen=G.GetMorganGenerator(radius=2,fpSize=2048)
def load(g):
    t=T[g][0][0]; df=pd.read_json(f'{RAW}/{t}.jsonl',lines=True)
    if df.empty: return df,{}
    df=df[(df.standard_relation=='=')&df.standard_value.notna()&df.canonical_smiles.notna()]
    df=df[df.data_validity_comment.isna()&(df.potential_duplicate.fillna(0).astype(int)==0)]
    df['pX']=9-np.log10(df.standard_value.astype(float).clip(lower=1e-6)); df=df[(df.pX>=2)&(df.pX<=12.5)]
    smi=df.drop_duplicates('molecule_chembl_id').set_index('molecule_chembl_id').canonical_smiles.to_dict()
    r=df.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().rename('pX').reset_index()
    return r,smi
fpc={}
def fp(m,smi):
    if m not in fpc:
        mol=Chem.MolFromSmiles(smi); fpc[m]=None if mol is None else gen.GetCountFingerprintAsNumPy(mol).astype(np.float32)
    return fpc[m]
def rf(Xtr,ytr,seed): return RandomForestRegressor(100,max_features='sqrt',min_samples_leaf=1,n_jobs=2,random_state=seed).fit(Xtr,ytr)
rows=[];anch=[]
for l in open('data/pairs.txt'):
    A,B=l.strip().split('/'); a,sa=load(A); b,sb=load(B)
    if a.empty or b.empty: continue
    smi={**sa,**sb}
    m=a.merge(b,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    if m.molecule_chembl_id.nunique()<60: continue
    res=dict(pair=f'{A}/{B}',n_rows=len(m),n_mols=m.molecule_chembl_id.nunique(),n_docs=m.document_chembl_id.nunique())
    # ---- scaffold split on per-compound median delta
    c=m.groupby('molecule_chembl_id').d.median().reset_index()
    def scaf(s):
        try: return MurckoScaffold.MurckoScaffoldSmiles(smi[s])
        except: return s
    c['scaf']=c.molecule_chembl_id.map(scaf)
    X=np.stack([fp(x,smi[x]) for x in c.molecule_chembl_id]); y=c.d.values
    scs=c.scaf.unique(); r2s=[];leak=[]
    docs_of=m.groupby('molecule_chembl_id').document_chembl_id.apply(set).to_dict()
    for seed in range(5):
        rng=np.random.default_rng(seed); te_s=set(rng.choice(scs,max(1,int(.25*len(scs))),replace=False))
        te=c.scaf.isin(te_s).values
        if te.sum()<10 or (~te).sum()<10: continue
        p=rf(X[~te],y[~te],seed).predict(X[te]); r2s.append(r2_score(y[te],p))
        trdocs=set().union(*[docs_of[x] for x in c.molecule_chembl_id[~te]])
        leak.append(np.mean([len(docs_of[x]&trdocs)>0 for x in c.molecule_chembl_id[te]]))
    res.update(r2_scaffold=np.mean(r2s) if r2s else np.nan,leak=np.mean(leak) if leak else np.nan)
    # ---- variance decomposition (docs with >=5 compounds)
    dc=m.groupby('document_chembl_id').size(); good=dc[dc>=5].index; mg=m[m.document_chembl_id.isin(good)].copy()
    if len(good)>=8 and mg.molecule_chembl_id.nunique()>=60:
        mu=mg.d.mean(); res['frac_between']=((mg.groupby('document_chembl_id').d.transform('mean')-mu)**2).sum()/((mg.d-mu)**2).sum()
        # ---- LPO
        Xr=np.stack([fp(x,smi[x]) for x in mg.molecule_chembl_id]); yr=mg.d.values; dd=mg.document_chembl_id.values
        mol=mg.molecule_chembl_id.values; udocs=np.array(sorted(set(dd)))
        R=dict(r2=[],r2w=[],rhow=[]); AN={k:[] for k in (0,1,2,3,5)}
        for seed in range(5):
            rng=np.random.default_rng(100+seed); te_d=set(rng.choice(udocs,max(1,int(.25*len(udocs))),replace=False))
            te=np.isin(dd,list(te_d)); temol=set(mol[te]); tr=~te&~np.isin(mol,list(temol))
            if te.sum()<20 or tr.sum()<30: continue
            p=rf(Xr[tr],yr[tr],seed).predict(Xr[te]); yt=yr[te]; dt=dd[te]
            R['r2'].append(r2_score(yt,p))
            pc=pd.Series(p).groupby(dt).transform(lambda s:s-s.mean()).values; yc=pd.Series(yt).groupby(dt).transform(lambda s:s-s.mean()).values
            R['r2w'].append(r2_score(yc,pc)); R['rhow'].append(spearmanr(yc,pc)[0])
            for k in AN:
                P=[];Y=[]
                for dname in set(dt):
                    ix=np.where(dt==dname)[0]
                    if k>0:
                        if len(ix)<k+3: continue
                        sel=rng.choice(ix,k,replace=False); off=np.mean(yt[sel]-p[sel]); rest=np.setdiff1d(ix,sel)
                    else: off=0; rest=ix
                    P.append(p[rest]+off); Y.append(yt[rest])
                if P: AN[k].append(r2_score(np.concatenate(Y),np.concatenate(P)))
        res.update(r2_lpo=np.mean(R['r2']) if R['r2'] else np.nan,r2_within=np.mean(R['r2w']) if R['r2w'] else np.nan,rho_within=np.nanmean(R['rhow']) if R['rhow'] else np.nan,
                   n_docs_lpo=len(good))
        for k in AN: res[f'r2_anchor{k}']=np.mean(AN[k]) if AN[k] else np.nan
    rows.append(res); print(res,flush=True)
    pd.DataFrame(rows).to_csv(f'{OUT}/models_by_pair.csv',index=False)
