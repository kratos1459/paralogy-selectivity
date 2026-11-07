"""Export the benchmark: per (pair, compound, document) selectivity values with publication-level (5 seeds) and chronological split labels,
plus compounds.csv.gz (SMILES). Run from repo root with raw/ present. Split rules: documents with >=5 matched compounds; pairs with >=8 such documents and >=60 compounds;
LPO: 25% of documents held out per seed, training = rows whose compound is absent from every test document ('excluded' otherwise);
time: documents after the 75th-percentile year are test, compounds appearing in test documents are excluded from training."""
import json,os,numpy as np,pandas as pd
T=json.load(open('data/targets_resolved.json')); yr=pd.read_csv('data/doc_years.csv').set_index('document_chembl_id').year.to_dict()
def load(g):
    t=T[g][0][0]; df=pd.read_json(f'raw/{t}.jsonl',lines=True)
    df=df[(df.standard_relation=='=')&df.standard_value.notna()&df.canonical_smiles.notna()&df.data_validity_comment.isna()&(df.potential_duplicate.fillna(0).astype(int)==0)]
    df['pX']=9-np.log10(df.standard_value.astype(float).clip(lower=1e-6)); df=df[(df.pX>=2)&(df.pX<=12.5)]
    return df.groupby(['molecule_chembl_id','document_chembl_id']).pX.median().rename('pX').reset_index(),df.drop_duplicates('molecule_chembl_id').set_index('molecule_chembl_id').canonical_smiles.to_dict()
out=[];smi={}
for l in open('data/pairs.txt'):
    A,B=l.strip().split('/'); a,sa=load(A); b,sb=load(B); smi.update(sa); smi.update(sb)
    m=a.merge(b,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['delta']=m.pXA-m.pXB
    dc=m.groupby('document_chembl_id').size(); m=m[m.document_chembl_id.isin(dc[dc>=5].index)].copy()
    if m.document_chembl_id.nunique()<8 or m.molecule_chembl_id.nunique()<60: continue
    dd=m.document_chembl_id.values; mol=m.molecule_chembl_id.values; ud=np.array(sorted(set(dd)))
    for s in range(5):
        rng=np.random.default_rng(100+s); td=set(rng.choice(ud,max(1,int(.25*len(ud))),replace=False)); te=np.isin(dd,list(td)); tm=set(mol[te])
        m[f'split_lpo_s{s}']=np.where(te,'test',np.where(np.isin(mol,list(tm)),'excluded','train'))
    m['year']=m.document_chembl_id.map(yr); cut=m.year.quantile(.75); te=(m.year>cut).values; tm=set(mol[te])
    m['split_time']=np.where(te,'test',np.where(np.isin(mol,list(tm)),'excluded','train')); m['time_cut_year']=cut
    m.insert(0,'pair',f'{A}/{B}'); out.append(m.rename(columns={'pXA':f'pX_target_A','pXB':'pX_target_B'}))
o=pd.concat(out); os.makedirs('benchmark',exist_ok=True)
o.round(4).to_csv('benchmark/selectivity_benchmark.csv.gz',index=False)
ids=sorted(set(o.molecule_chembl_id)); pd.DataFrame({'molecule_chembl_id':ids,'canonical_smiles':[smi[i] for i in ids]}).to_csv('benchmark/compounds.csv.gz',index=False)
print(len(o),'rows',o.pair.nunique(),'pairs',len(ids),'compounds')
