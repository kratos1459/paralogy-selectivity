"""XGBoost (300 trees, depth 6) versus random forest under leave-publication-out, same capped training as 09."""
import sys,json,os,numpy as np,pandas as pd
from rdkit import Chem,RDLogger; from rdkit.Chem import rdFingerprintGenerator as G
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor,HistGradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import GroupKFold
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
    df['cls']=np.where(df.standard_type=='IC50','IC50','KiKd')
    r=df.groupby(['molecule_chembl_id','document_chembl_id','cls']).pX.median().rename('pX').reset_index()
    return r,smi
fpc={}
def fp(m,smi):
    if m not in fpc:
        mol=Chem.MolFromSmiles(smi); fpc[m]=None if mol is None else gen.GetCountFingerprintAsNumPy(mol).astype(np.float32)
    return fpc[m]
def rf(Xtr,ytr,seed): return RandomForestRegressor(100,max_features='sqrt',min_samples_leaf=1,n_jobs=2,random_state=seed).fit(Xtr,ytr)

import xgboost as xgb
CAP=2500
def within_r2(yt,p,dt):
    pc=pd.Series(p).groupby(dt).transform(lambda s:s-s.mean()).values; yc=pd.Series(yt).groupby(dt).transform(lambda s:s-s.mean()).values
    return r2_score(yc,pc)
rows=[]
for l in open('data/pairs.txt'):
    A,B=l.strip().split('/'); a,sa=load(A); b,sb=load(B)
    if a.empty or b.empty: continue
    smi={**sa,**sb}
    a2=a.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().reset_index(); b2=b.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().reset_index()
    m=a2.merge(b2,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    if m.empty: continue
    dc=m.groupby('document_chembl_id').size(); good=dc[dc>=5].index; mg=m[m.document_chembl_id.isin(good)].copy()
    if not (len(good)>=8 and mg.molecule_chembl_id.nunique()>=60): continue
    Xr=np.stack([fp(x,smi[x]) for x in mg.molecule_chembl_id]); yr=mg.d.values; dd=mg.document_chembl_id.values
    mol=mg.molecule_chembl_id.values; udocs=np.array(sorted(set(dd)))
    acc={k:[] for k in ('rf_p','rf_w','xgb_p','xgb_w','ens_p','ens_w')}
    for seed in range(3):
        rng=np.random.default_rng(500+seed); te_d=set(rng.choice(udocs,max(1,int(.25*len(udocs))),replace=False))
        te=np.isin(dd,list(te_d)); temol=set(mol[te]); tr=~te&~np.isin(mol,list(temol))
        if te.sum()<20 or tr.sum()<30: continue
        ti=np.where(tr)[0]
        if len(ti)>CAP: ti=rng.choice(ti,CAP,replace=False)
        Xt,yt_=Xr[ti],yr[ti]; Xe,ye,de=Xr[te],yr[te],dd[te]
        p1=rf(Xt,yt_,seed).predict(Xe)
        p2=xgb.XGBRegressor(n_estimators=300,max_depth=6,learning_rate=0.1,subsample=0.8,colsample_bytree=0.5,n_jobs=2,random_state=seed,verbosity=0).fit(Xt,yt_).predict(Xe)
        for k,p in (('rf',p1),('xgb',p2),('ens',(p1+p2)/2)): acc[k+'_p'].append(r2_score(ye,p)); acc[k+'_w'].append(within_r2(ye,p,de))
    r=dict(pair=f'{A}/{B}',**{k:(np.mean(v) if v else np.nan) for k,v in acc.items()})
    rows.append(r); print(r['pair'],round(r['rf_p'],2),round(r['xgb_p'],2),flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/xgboost_vs_rf_by_pair.csv',index=False)
