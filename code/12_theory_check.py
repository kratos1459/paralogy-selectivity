"""Check the closed-form shrinkage theory: with offset variance tau2 across publications and residual variance sigma2, the plain offset has excess error sigma2/m and the empirical-Bayes offset w*sigma2/m with w = tau2/(tau2+sigma2/m). Compares predicted with observed MSE per pair."""
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

CAP=6000
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
    acc={f'{s_}{k}':[] for s_ in ('off','eb') for k in (1,2,3,5)}; acc['r0']=[]
    for k in (1,2,3,5):
        for nm in ('mse_off','mse_eb','mse_off_th','mse_eb_th'): acc[f'{nm}{k}']=[]
    for nm in ('tau2','s2tr','sig2te'): acc[nm]=[]
    for seed in range(3):
        rng=np.random.default_rng(400+seed); te_d=set(rng.choice(udocs,max(1,int(.25*len(udocs))),replace=False))
        te=np.isin(dd,list(te_d)); temol=set(mol[te]); tr=~te&~np.isin(mol,list(temol))
        if te.sum()<20 or tr.sum()<30: continue
        ti=np.where(tr)[0]
        if len(ti)>CAP: ti=rng.choice(ti,CAP,replace=False)
        Xt,yt_,dtr=Xr[ti],yr[ti],dd[ti]; Xe,ye,de=Xr[te],yr[te],dd[te]
        p=RandomForestRegressor(60,max_features='sqrt',n_jobs=2,random_state=seed).fit(Xt,yt_).predict(Xe); acc['r0'].append(r2_score(ye,p))
        if len(np.unique(dtr))<6: continue
        oof=np.full(len(ti),np.nan); mt=mol[ti]
        for a_,b_ in GroupKFold(3).split(Xt,yt_,dtr):
            bm=set(mt[b_]); keep=np.array([x not in bm for x in mt[a_]])
            oof[b_]=RandomForestRegressor(60,max_features='sqrt',n_jobs=2,random_state=seed).fit(Xt[a_[keep]],yt_[a_[keep]]).predict(Xt[b_])
        res=yt_-oof; g=pd.DataFrame({'d':dtr,'r':res}).dropna().groupby('d').r.agg(['mean','var','size']); g=g[g['size']>=2]
        s2=np.average(g['var'],weights=g['size']-1); mu=np.average(g['mean'],weights=g['size'])
        tau2=max(np.average((g['mean']-mu)**2,weights=g['size'])-np.mean(s2/g['size']),0.0)
        rte=ye-p; sig2te=np.average([np.var(rte[de==u]) for u in set(de)],weights=[np.sum(de==u) for u in set(de)])
        acc['tau2'].append(tau2); acc['s2tr'].append(s2); acc['sig2te'].append(sig2te)
        for k in (1,2,3,5):
            Po={'off':[],'eb':[]};Y=[]
            for dn in set(de):
                ix=np.where(de==dn)[0]
                if len(ix)<k+3: continue
                sel=rng.choice(ix,k,replace=False); rest=np.setdiff1d(ix,sel); raw=np.mean(ye[sel]-p[sel]); w=tau2/(tau2+s2/k)
                Po['off'].append(p[rest]+raw); Po['eb'].append(p[rest]+mu+w*(raw-mu)); Y.append(ye[rest])
            if Y:
                Yc=np.concatenate(Y)
                for nm in Po:
                    Pc=np.concatenate(Po[nm]); acc[f'{nm}{k}'].append(r2_score(Yc,Pc)); acc['mse_'+nm+str(k)].append(np.mean((Yc-Pc)**2))
                w_=tau2/(tau2+s2/k); acc[f'mse_off_th{k}'].append(sig2te*(1+1/k)); acc[f'mse_eb_th{k}'].append(sig2te*(1+w_/k))
    r=dict(pair=f'{A}/{B}',family=None,**{k:(np.mean(v) if v else np.nan) for k,v in acc.items()})
    rows.append(r); print(r['pair'],flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/theory_check_by_pair.csv',index=False)
