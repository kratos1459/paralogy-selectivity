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

rows=[]
for l in open('data/pairs.txt'):
    A,B=l.strip().split('/'); a,sa=load(A); b,sb=load(B)
    if a.empty or b.empty: continue
    smi={**sa,**sb}
    m=a.merge(b,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    if m.molecule_chembl_id.nunique()<60: continue
    dc=m.groupby('document_chembl_id').size(); good=dc[dc>=5].index; mg=m[m.document_chembl_id.isin(good)].copy()
    if not (len(good)>=8 and mg.molecule_chembl_id.nunique()>=60): continue
    Xr=np.stack([fp(x,smi[x]) for x in mg.molecule_chembl_id]); yr=mg.d.values; dd=mg.document_chembl_id.values
    mol=mg.molecule_chembl_id.values; udocs=np.array(sorted(set(dd)))
    out={(k,s):[] for k in (1,2,3,5) for s in ('offset','shrunk')}
    for seed in range(5):
        rng=np.random.default_rng(100+seed); te_d=set(rng.choice(udocs,max(1,int(.25*len(udocs))),replace=False))
        te=np.isin(dd,list(te_d)); temol=set(mol[te]); tr=~te&~np.isin(mol,list(temol))
        if te.sum()<20 or tr.sum()<30: continue
        mdl=RandomForestRegressor(100,max_features='sqrt',n_jobs=2,random_state=seed,oob_score=True).fit(Xr[tr],yr[tr])
        p=mdl.predict(Xr[te]); yt=yr[te]; dt=dd[te]
        res=yr[tr]-mdl.oob_prediction_; dtr=dd[tr]; ok=np.isfinite(res)
        g=pd.DataFrame({'d':dtr[ok],'r':res[ok]}).groupby('d').r.agg(['mean','var','size']); g=g[g['size']>=2]
        s2=np.average(g['var'],weights=g['size']-1); mu=np.average(g['mean'],weights=g['size'])
        tau2=max(np.average((g['mean']-mu)**2,weights=g['size'])-np.mean(s2/g['size']),0.0)
        for k in (1,2,3,5):
            P={'offset':[],'shrunk':[]};Y=[]
            for dname in set(dt):
                ix=np.where(dt==dname)[0]
                if len(ix)<k+3: continue
                sel=rng.choice(ix,k,replace=False); rest=np.setdiff1d(ix,sel); raw=np.mean(yt[sel]-p[sel])
                w=tau2/(tau2+s2/k) if tau2+s2>0 else 0
                P['offset'].append(p[rest]+raw); P['shrunk'].append(p[rest]+mu+w*(raw-mu)); Y.append(yt[rest])
            if Y:
                Yc=np.concatenate(Y)
                for sname in P: out[(k,sname)].append(r2_score(Yc,np.concatenate(P[sname])))
    r=dict(pair=f'{A}/{B}')
    for (k,sname),v in out.items(): r[f'r2_{sname}{k}']=np.mean(v) if v else np.nan
    rows.append(r); print(r,flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/anchor_shrinkage_by_pair.csv',index=False)
