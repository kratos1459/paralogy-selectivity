"""Separate measurable protocol factors from residual disagreement using ChEMBL assay descriptors (run from repo root; needs raw/ and assays.csv from 13_fetch_assays step in README).
For compounds measured on both targets in >=2 documents, each document pair gives a squared difference of delta; half of its mean is Vx.
Document pairs are 'concordant' if, for both targets, the modal assay descriptor (standard type, BAO format, assay type, cell type present/absent, organism) matches;
otherwise 'discordant'. The share of Vx that disappears under concordance estimates protocol variance explainable by recorded descriptors."""
import json,itertools,numpy as np,pandas as pd
T=json.load(open('data/targets_resolved.json')); A=pd.read_csv('assays.csv').drop_duplicates('assay_chembl_id').set_index('assay_chembl_id')
A['cell']=A.assay_cell_type.notna().astype(int)
A['desc']=A.assay_type.astype(str)+'|'+A.bao_format.astype(str)+'|'+A.cell.astype(str)+'|'+A.assay_organism.astype(str)
def load(g):
    t=T[g][0][0]; df=pd.read_json(f'raw/{t}.jsonl',lines=True)
    df=df[(df.standard_relation=='=')&df.standard_value.notna()&df.data_validity_comment.isna()&(df.potential_duplicate.fillna(0).astype(int)==0)]
    df['pX']=9-np.log10(df.standard_value.astype(float).clip(lower=1e-6)); df=df[(df.pX>=2)&(df.pX<=12.5)]
    df['desc']=df.standard_type.astype(str)+'|'+df.assay_chembl_id.map(A.desc).fillna('NA')
    k=['molecule_chembl_id','document_chembl_id']
    pX=df.groupby(k).pX.median(); d=df.groupby(k).desc.agg(lambda s:s.value_counts().index[0])
    return pd.concat([pX,d],axis=1).reset_index()
rows=[]
for l in open('data/pairs.txt'):
    a,b=l.strip().split('/'); X=load(a); Y=load(b)
    m=X.merge(Y,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B')); m['d']=m.pXA-m.pXB
    m=m.groupby('molecule_chembl_id').filter(lambda x:x.document_chembl_id.nunique()>=2)
    conc=[];disc=[]
    for _,g in m.groupby('molecule_chembl_id'):
        r=g.to_dict('records')
        for u,v in itertools.combinations(r,2):
            sq=(u['d']-v['d'])**2/2
            (conc if (u['descA']==v['descA'] and u['descB']==v['descB']) else disc).append(sq)
    if len(conc)>=10 and len(disc)>=10: rows.append(dict(pair=f'{a}/{b}',n_conc=len(conc),n_disc=len(disc),Vx_conc=np.mean(conc),Vx_disc=np.mean(disc)))
    print(a,b,len(conc),len(disc),flush=True)
o=pd.DataFrame(rows); o['explained']=1-o.Vx_conc/o.Vx_disc; o.to_csv('results/assay_descriptor_concordance.csv',index=False)
w=lambda c:np.average(o[c],weights=o.n_conc+o.n_disc)
print(len(o),'pairs; pooled Vx concordant',round(w('Vx_conc'),3),'discordant',round(w('Vx_disc'),3),'median explained',round(o.explained.median(),3),'pairs with concordant<discordant',(o.Vx_conc<o.Vx_disc).sum())
