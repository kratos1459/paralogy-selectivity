"""External-style validation: hold out each large single-source document (>=150 matched compounds on a pair; typically a profiling panel from one source) and evaluate the model, plain offset and empirical-Bayes offset."""
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

CAP=6000; MINROWS=150
rows=[]
for l in open('data/pairs.txt'):
    A,B=l.strip().split('/'); a,sa=load(A); b,sb=load(B)
    if a.empty or b.empty: continue
    smi={**sa,**sb}
    a2=a.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().reset_index(); b2=b.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().reset_index()
    m=a2.merge(b2,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    if m.empty: continue
    sizes=m.groupby('document_chembl_id').size(); big=sizes[sizes>=MINROWS].index
    for dn in big:
        te=(m.document_chembl_id==dn).values; temol=set(m.molecule_chembl_id[te]); tr=(~te)&(~m.molecule_chembl_id.isin(temol)).values
        if tr.sum()<100: continue
        X=np.stack([fp(x,smi[x]) for x in m.molecule_chembl_id]); y=m.d.values; dd=m.document_chembl_id.values; mol=m.molecule_chembl_id.values
        rng=np.random.default_rng(7); ti=np.where(tr)[0]
        if len(ti)>CAP: ti=rng.choice(ti,CAP,replace=False)
        Xt,yt_,dtr=X[ti],y[ti],dd[ti]; Xe,ye=X[te],y[te]
        p=RandomForestRegressor(60,max_features='sqrt',n_jobs=2,random_state=0).fit(Xt,yt_).predict(Xe)
        r=dict(pair=f'{A}/{B}',doc=dn,n_test=int(te.sum()),n_train=len(ti),r2=r2_score(ye,p),r2_within=r2_score(ye-ye.mean(),p-p.mean()),rho=spearmanr(ye,p)[0],ceiling_note='single large document')
        if len(np.unique(dtr))>=6:
            oof=np.full(len(ti),np.nan); mt=mol[ti]
            for a_,b_ in GroupKFold(3).split(Xt,yt_,dtr):
                bm=set(mt[b_]); keep=np.array([x not in bm for x in mt[a_]])
                oof[b_]=RandomForestRegressor(60,max_features='sqrt',n_jobs=2,random_state=0).fit(Xt[a_[keep]],yt_[a_[keep]]).predict(Xt[b_])
            res=yt_-oof; g=pd.DataFrame({'d':dtr,'r':res}).dropna().groupby('d').r.agg(['mean','var','size']); g=g[g['size']>=2]
            s2=np.average(g['var'],weights=g['size']-1); mu=np.average(g['mean'],weights=g['size'])
            tau2=max(np.average((g['mean']-mu)**2,weights=g['size'])-np.mean(s2/g['size']),0.0)
            for k in (1,2,3,5):
                o1=[];o2=[]
                for rep in range(30):
                    sel=rng.choice(len(ye),k,replace=False); rest=np.setdiff1d(np.arange(len(ye)),sel); raw=np.mean(ye[sel]-p[sel]); w=tau2/(tau2+s2/k)
                    o1.append(r2_score(ye[rest],p[rest]+raw)); o2.append(r2_score(ye[rest],p[rest]+mu+w*(raw-mu)))
                r[f'off{k}']=np.mean(o1); r[f'eb{k}']=np.mean(o2)
        rows.append(r); print(r['pair'],dn,r['n_test'],round(r['r2'],2),flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/single_document_holdout.csv',index=False)
