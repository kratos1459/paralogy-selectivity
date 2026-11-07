"""Direct selectivity RF vs difference of two single-target RF potency models (scaffold split, 5 seeds). python 22 RAW OUT"""
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
    if m.empty: continue
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    if m.molecule_chembl_id.nunique()<60: continue
    c=m.groupby('molecule_chembl_id').agg(d=('d','median'),a=('pXA','median'),b=('pXB','median')).reset_index()
    def scaf(s_):
        try: return MurckoScaffold.MurckoScaffoldSmiles(smi[s_])
        except: return s_
    c['scaf']=c.molecule_chembl_id.map(scaf); X=np.stack([fp(x,smi[x]) for x in c.molecule_chembl_id]); scs=c.scaf.unique(); rd=[];rs=[]
    for seed in range(5):
        rng=np.random.default_rng(seed); te_s=set(rng.choice(scs,max(1,int(.25*len(scs))),replace=False)); te=c.scaf.isin(te_s).values
        if te.sum()<10 or (~te).sum()<10: continue
        y=c.d.values; pd_=rf(X[~te],y[~te],seed).predict(X[te])
        pa=rf(X[~te],c.a.values[~te],seed).predict(X[te]); pb=rf(X[~te],c.b.values[~te],seed).predict(X[te])
        rd.append(r2_score(y[te],pd_)); rs.append(r2_score(y[te],pa-pb))
    if rd: rows.append(dict(pair=f'{A}/{B}',r2_direct=np.mean(rd),r2_diff_of_two=np.mean(rs))); print(rows[-1],flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/direct_vs_difference_by_pair.csv',index=False)
