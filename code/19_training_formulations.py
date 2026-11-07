"""Standard vs series-centred vs intra-publication pairwise training x RF/XGBoost/MLP on leave-publication-out splits.
Outputs within-series R2, Spearman rho, dispersion ratio k, centred Pearson r, pooled R2. Usage: python 19 RAW OUT [part nparts]"""
import sys,json,os,numpy as np,pandas as pd,xgboost as xgb
from rdkit import Chem,RDLogger; from rdkit.Chem import rdFingerprintGenerator as G
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import r2_score
from scipy.stats import spearmanr,pearsonr
RDLogger.DisableLog('rdApp.*')
RAW=sys.argv[1]; OUT=sys.argv[2]; part=int(sys.argv[3]) if len(sys.argv)>3 else 0; nparts=int(sys.argv[4]) if len(sys.argv)>4 else 1
os.makedirs(OUT,exist_ok=True)
T=json.load(open('data/targets_resolved.json')); T['ADRA1A']=[['CHEMBL229','a']]
gen=G.GetMorganGenerator(radius=2,fpSize=2048)
def load(g):
    t=T[g][0][0]; df=pd.read_json(f'{RAW}/{t}.jsonl',lines=True)
    if df.empty: return df,{}
    df=df[(df.standard_relation=='=')&df.standard_value.notna()&df.canonical_smiles.notna()]
    df=df[df.data_validity_comment.isna()&(df.potential_duplicate.fillna(0).astype(int)==0)]
    df['pX']=9-np.log10(df.standard_value.astype(float).clip(lower=1e-6)); df=df[(df.pX>=2)&(df.pX<=12.5)]
    smi=df.drop_duplicates('molecule_chembl_id').set_index('molecule_chembl_id').canonical_smiles.to_dict()
    return df.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().reset_index(),smi
fpc={}
def fp(m,smi):
    if m not in fpc:
        mol=Chem.MolFromSmiles(smi); fpc[m]=None if mol is None else gen.GetCountFingerprintAsNumPy(mol).astype(np.float32)
    return fpc[m]
def mk(name,seed):
    if name=='rf': return RandomForestRegressor(100,max_features='sqrt',n_jobs=1,random_state=seed)
    if name=='xgb': return xgb.XGBRegressor(n_estimators=300,max_depth=6,learning_rate=0.1,subsample=0.8,colsample_bytree=0.5,n_jobs=1,random_state=seed,verbosity=0)
    return MLPRegressor(hidden_layer_sizes=(256,64),max_iter=40,random_state=seed)
def tf(name,X): return X if name!='mlp' else np.log1p(np.abs(X))*np.sign(X)
def centre(v,d): return v-pd.Series(v).groupby(d).transform('mean').values
CAP=2000; NPAIR=4000; NREF=20
rows=[]
pairs=[l.strip() for l in open('data/pairs.txt')]
for pi,l in enumerate(pairs):
    if pi%nparts!=part: continue
    A,B=l.split('/'); a,sa=load(A); b,sb=load(B)
    if a.empty or b.empty: continue
    smi={**sa,**sb}; m=a.merge(b,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    if m.empty: continue
    dc=m.groupby('document_chembl_id').size(); good=dc[dc>=5].index; mg=m[m.document_chembl_id.isin(good)].copy()
    if not (len(good)>=8 and mg.molecule_chembl_id.nunique()>=60): continue
    Xr=np.stack([fp(x,smi[x]) for x in mg.molecule_chembl_id]); yr=mg.d.values; dd=mg.document_chembl_id.values; mol=mg.molecule_chembl_id.values; ud=np.array(sorted(set(dd)))
    acc={}
    for seed in range(3):
        rng=np.random.default_rng(200+seed); te_d=set(rng.choice(ud,max(1,int(.25*len(ud))),replace=False))
        te=np.isin(dd,list(te_d)); tr=~te&~np.isin(mol,list(set(mol[te])))
        if te.sum()<20 or tr.sum()<30: continue
        ti=np.where(tr)[0]
        if len(ti)>CAP: ti=rng.choice(ti,CAP,replace=False)
        Xt,yt,dt=Xr[ti],yr[ti],dd[ti]; Xe,ye,de=Xr[te],yr[te],dd[te]; yec=centre(ye,de)
        # candidate within-document training pairs
        idx=[np.where(dt==u)[0] for u in np.unique(dt)]; idx=[i for i in idx if len(i)>=2]
        pr=[]
        for _ in range(NPAIR):
            i=idx[rng.integers(len(idx))]; j=rng.choice(i,2,replace=False); pr.append(j)
        pr=np.array(pr)
        for name in ('rf','xgb','mlp'):
            def rec(scheme,pred_c,pooled=None):
                k=np.std(pred_c)/np.std(yec); r=pearsonr(pred_c,yec)[0] if np.std(pred_c)>0 else 0.0
                d=acc.setdefault((scheme,name),{'w':[],'rho':[],'k':[],'r':[],'p':[]})
                d['w'].append(r2_score(yec,pred_c)); d['rho'].append(spearmanr(pred_c,yec)[0]); d['k'].append(k); d['r'].append(r)
                if pooled is not None: d['p'].append(r2_score(ye,pooled))
            p=mk(name,seed).fit(tf(name,Xt),yt).predict(tf(name,Xe)); rec('standard',centre(p,de),p)
            ytc=centre(yt,dt); p=mk(name,seed).fit(tf(name,Xt),ytc).predict(tf(name,Xe))+yt.mean(); rec('series_centred',centre(p,de),p)
            mdl=mk(name,seed).fit(tf(name,Xt[pr[:,0]]-Xt[pr[:,1]]),yt[pr[:,0]]-yt[pr[:,1]]); sc=np.zeros(len(ye))
            for u in np.unique(de):
                ix=np.where(de==u)[0]
                for i in ix:
                    oth=rng.choice(ix[ix!=i],min(NREF,len(ix)-1),replace=False)
                    sc[i]=mdl.predict(tf(name,Xe[i]-Xe[oth])).mean()
            rec('pairwise',centre(sc,de))
    for (scheme,name),d in acc.items():
        rows.append(dict(pair=f'{A}/{B}',scheme=scheme,model=name,**{k:np.mean(v) if v else np.nan for k,v in d.items()}))
    print(A,B,len(rows),flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/training_formulations_part{part}.csv',index=False)
