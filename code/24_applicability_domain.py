"""Applicability domain: max Tanimoto similarity of LPO test compounds to the training set versus absolute error. python 24 RAW OUT"""
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
from rdkit import DataStructs
rows=[]
for l in open('data/pairs.txt'):
    A,B=l.strip().split('/'); a,sa=load(A); b,sb=load(B)
    if a.empty or b.empty: continue
    smi={**sa,**sb}; m=a.merge(b,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    if m.empty: continue
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    dc=m.groupby('document_chembl_id').size(); good=dc[dc>=5].index; mg=m[m.document_chembl_id.isin(good)].copy()
    if not (len(good)>=8 and mg.molecule_chembl_id.nunique()>=60): continue
    Xr=np.stack([fp(x,smi[x]) for x in mg.molecule_chembl_id]); yr=mg.d.values; dd=mg.document_chembl_id.values; mol=mg.molecule_chembl_id.values; udocs=np.array(sorted(set(dd)))
    bits={x:DataStructs.CreateFromBitString(''.join('1' if v>0 else '0' for v in fp(x,smi[x]))) for x in set(mol)}
    for seed in range(2):
        rng=np.random.default_rng(100+seed); te_d=set(rng.choice(udocs,max(1,int(.25*len(udocs))),replace=False))
        te=np.isin(dd,list(te_d)); temol=set(mol[te]); tr=~te&~np.isin(mol,list(temol))
        if te.sum()<20 or tr.sum()<30: continue
        ti=np.where(tr)[0]
        if len(ti)>3000: ti=rng.choice(ti,3000,replace=False)
        p=rf(Xr[ti],yr[ti],seed).predict(Xr[te]); trb=[bits[x] for x in sorted(set(mol[ti]))]
        sim=np.array([max(DataStructs.BulkTanimotoSimilarity(bits[x],trb)) for x in mol[te]])
        yc=yr[te]-pd.Series(yr[te]).groupby(dd[te]).transform('mean').values; pc=p-pd.Series(p).groupby(dd[te]).transform('mean').values
        for s_,e_,ec in zip(sim,np.abs(yr[te]-p),np.abs(yc-pc)): rows.append(dict(pair=f'{A}/{B}',seed=seed,maxsim=s_,abs_err=e_,abs_err_centred=ec))
    print(A,B,len(rows),flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/applicability_domain_records.csv',index=False)
