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

CAP=2500
def within_r2(yt,p,dt):
    pc=pd.Series(p).groupby(dt).transform(lambda s:s-s.mean()).values; yc=pd.Series(yt).groupby(dt).transform(lambda s:s-s.mean()).values
    return r2_score(yc,pc)
rows=[];strat=[]
for l in open('data/pairs.txt'):
    A,B=l.strip().split('/'); a,sa=load(A); b,sb=load(B)
    if a.empty or b.empty: continue
    smi={**sa,**sb}
    # ---- stratified between-publication share (equilibrium Ki/Kd vs functional IC50), same class for both targets
    for cls in ('KiKd','IC50'):
        if a[a.cls==cls].empty or b[b.cls==cls].empty: continue
        mm=a[a.cls==cls].merge(b[b.cls==cls],on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); mm['d']=mm.pXA-mm.pXB
        dc=mm.groupby('document_chembl_id').size(); mm=mm[mm.document_chembl_id.isin(dc[dc>=5].index)]
        if mm.document_chembl_id.nunique()>=5 and len(mm)>=60:
            mu=mm.d.mean(); strat.append(dict(pair=f'{A}/{B}',cls=cls,n=len(mm),docs=mm.document_chembl_id.nunique(),frac_between=((mm.groupby('document_chembl_id').d.transform('mean')-mu)**2).sum()/((mm.d-mu)**2).sum()))
    pd.DataFrame(strat).to_csv(f'{OUT}/stratified_variance_share.csv',index=False)
    a2=a.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().reset_index(); b2=b.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().reset_index()
    m=a2.merge(b2,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    if m.empty: continue
    dc=m.groupby('document_chembl_id').size(); good=dc[dc>=5].index; mg=m[m.document_chembl_id.isin(good)].copy()
    if not (len(good)>=8 and mg.molecule_chembl_id.nunique()>=60): continue
    Xr=np.stack([fp(x,smi[x]) for x in mg.molecule_chembl_id]); yr=mg.d.values; dd=mg.document_chembl_id.values
    mol=mg.molecule_chembl_id.values; udocs=np.array(sorted(set(dd)))
    acc={k:[] for k in ['rf_p','rf_w','hgb_p','hgb_w','mlp_p','mlp_w','cen_p','cen_w','ens_p','ens_w']+[f'{s}{k}' for s in ('off','eb') for k in (1,2,3,5)]}
    for seed in range(3):
        rng=np.random.default_rng(200+seed); te_d=set(rng.choice(udocs,max(1,int(.25*len(udocs))),replace=False))
        te=np.isin(dd,list(te_d)); temol=set(mol[te]); tr=~te&~np.isin(mol,list(temol))
        if te.sum()<20 or tr.sum()<30: continue
        ti=np.where(tr)[0]
        if len(ti)>CAP: ti=rng.choice(ti,CAP,replace=False)
        Xt,yt_,dtr=Xr[ti],yr[ti],dd[ti]; Xe,ye,de=Xr[te],yr[te],dd[te]
        P={}
        P['rf']=rf(Xt,yt_,seed).predict(Xe)
        P['hgb']=HistGradientBoostingRegressor(max_iter=200,learning_rate=0.1,random_state=seed).fit(np.log1p(Xt),yt_).predict(np.log1p(Xe))
        P['mlp']=MLPRegressor(hidden_layer_sizes=(256,64),max_iter=40,early_stopping=False,random_state=seed).fit(np.log1p(Xt),yt_).predict(np.log1p(Xe))
        cen=yt_-pd.Series(yt_).groupby(dtr).transform('mean').values
        P['cen']=rf(Xt,cen,seed).predict(Xe)+yt_.mean()
        P['ens']=(P['rf']+P['hgb']+P['mlp'])/3
        for k,p in P.items(): acc[k+'_p'].append(r2_score(ye,p)); acc[k+'_w'].append(within_r2(ye,p,de)) if k!='cen' else acc['cen_w'].append(within_r2(ye,p,de))
        # empirical Bayes anchor offset with publication-level nested CV for tau^2
        ug=np.unique(dtr)
        if len(ug)>=6:
            oof=np.full(len(ti),np.nan)
            for a_,b_ in GroupKFold(3).split(Xt,yt_,dtr):
                a_m=set(mol[ti[a_]]) ; keep=np.array([mol[ti[i]] not in set(mol[ti[b_]]) for i in a_])
                oof[b_]=rf(Xt[a_[keep]],yt_[a_[keep]],seed).predict(Xt[b_])
            res=yt_-oof; g=pd.DataFrame({'d':dtr,'r':res}).dropna().groupby('d').r.agg(['mean','var','size']); g=g[g['size']>=2]
            s2=np.average(g['var'],weights=g['size']-1); mu=np.average(g['mean'],weights=g['size'])
            tau2=max(np.average((g['mean']-mu)**2,weights=g['size'])-np.mean(s2/g['size']),0.0)
            for k in (1,2,3,5):
                Po={'off':[],'eb':[]};Y=[]
                for dn in set(de):
                    ix=np.where(de==dn)[0]
                    if len(ix)<k+3: continue
                    sel=rng.choice(ix,k,replace=False); rest=np.setdiff1d(ix,sel); raw=np.mean(ye[sel]-P['rf'][sel]); w=tau2/(tau2+s2/k)
                    Po['off'].append(P['rf'][rest]+raw); Po['eb'].append(P['rf'][rest]+mu+w*(raw-mu)); Y.append(ye[rest])
                if Y:
                    Yc=np.concatenate(Y)
                    for nm in Po: acc[f'{nm}{k}'].append(r2_score(Yc,np.concatenate(Po[nm])))
    r=dict(pair=f'{A}/{B}',**{k:(np.mean(v) if v else np.nan) for k,v in acc.items()})
    rows.append(r); print(r['pair'],round(r['rf_p'],2),round(r['hgb_p'],2),round(r['mlp_p'],2),flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/extra_models_by_pair.csv',index=False)
