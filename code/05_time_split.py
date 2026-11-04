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

yr=pd.read_csv('data/doc_years.csv').set_index('document_chembl_id').year.to_dict()
rows=[]
for l in open('data/pairs.txt'):
    A,B=l.strip().split('/'); a,sa=load(A); b,sb=load(B)
    if a.empty or b.empty: continue
    smi={**sa,**sb}
    m=a.merge(b,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    m=m[m.molecule_chembl_id.map(lambda x: fp(x,smi[x]) is not None)]
    if m.empty: continue
    dc=m.groupby('document_chembl_id').size(); m=m[m.document_chembl_id.isin(dc[dc>=5].index)].copy()
    if m.empty: continue
    m['year']=m.document_chembl_id.map(yr); m=m.dropna(subset=['year'])
    if m.document_chembl_id.nunique()<8 or m.molecule_chembl_id.nunique()<60: continue
    cut=m.year.quantile(.75); te=(m.year>cut).values
    # keep strictly later docs; need enough
    mol=m.molecule_chembl_id.values; tr=~te&~np.isin(mol,list(set(mol[te])))
    if te.sum()<30 or tr.sum()<60: continue
    X=np.stack([fp(x,smi[x]) for x in m.molecule_chembl_id]); y=m.d.values; dd=m.document_chembl_id.values
    p=rf(X[tr],y[tr],0).predict(X[te]); yt=y[te]; dt=dd[te]
    pc=pd.Series(p).groupby(dt).transform(lambda s:s-s.mean()).values; yc=pd.Series(yt).groupby(dt).transform(lambda s:s-s.mean()).values
    r=dict(pair=f'{A}/{B}',cut_year=cut,n_train=int(tr.sum()),n_test=int(te.sum()),r2_time=r2_score(yt,p),r2_within_time=r2_score(yc,pc),rho_within_time=spearmanr(yc,pc)[0])
    rows.append(r); print(r,flush=True); pd.DataFrame(rows).to_csv(f'{OUT}/time_split_by_pair.csv',index=False)
