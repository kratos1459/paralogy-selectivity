"""Variance-component analysis of paralog selectivity (delta = pX(A) - pX(B)) across ChEMBL publications.
Reads raw/<target>.jsonl (from 01_fetch_chembl.py) and data/targets_resolved.json.

Curation: relation '=', no data_validity_comment, potential_duplicate==0, 2 <= pX <= 12.5, pX = 9 - log10(nM);
median pX per (compound, document, target). delta is formed only within the same document.

Quantities per pair
  Vtot   total variance of delta over (compound, document) rows
  Vbet   variance of document means (between-publication component); share = Vbet/Vtot
  Vx     cross-document replicate variance: within-compound variance of delta across >=2 documents
         (compound-specific selectivity cancels, so Vx = protocol offset + measurement noise per measurement)
  ceiling = 1 - Vx/Vtot : upper bound on R^2 for predicting a held-out publication's delta
  prot_share_of_between = min(Vx/Vbet,1);  series_resid = max(Vbet-Vx,0)/Vtot
Bootstrap (500x, resampling compounds) gives 95% CIs for Vx and ceiling. Pairs need >=50 same-document rows;
Vx is reported only with >=15 repeat compounds, and summaries use pairs with >=20.
Caveat: repeat compounds are mostly reference/tool compounds, so Vx may not represent the typical series compound."""
import json,os,numpy as np,pandas as pd
rng=np.random.default_rng(0)
T=json.load(open('data/targets_resolved.json'))
def load(t):
    fn=f'raw/{t}.jsonl'
    df=pd.read_json(fn,lines=True)
    if df.empty: return df
    df=df[(df.standard_relation=='=')&df.standard_value.notna()&df.canonical_smiles.notna()]
    df=df[df.data_validity_comment.isna()&(df.potential_duplicate.fillna(0).astype(int)==0)]
    df['pX']=9-np.log10(df.standard_value.astype(float).clip(lower=1e-6))
    return df[(df.pX>=2)&(df.pX<=12.5)]
def curate(df):
    g=df.groupby(['molecule_chembl_id','document_chembl_id'])
    r=g.pX.median().rename('pX').reset_index(); r['nrec']=g.pX.size().values
    return r
def vcross(c):
    g=c.groupby('molecule_chembl_id').d
    ss=((c.d-g.transform('mean'))**2).sum()
    df=(g.size()-1).sum()
    return ss/df
def pair_stats(a,b,nboot=500):
    A,B=curate(a),curate(b)
    m=A.merge(B,on=['molecule_chembl_id','document_chembl_id'],suffixes=('A','B'))
    if len(m)<50: return None
    m['d']=m.pXA-m.pXB; mu=m.d.mean(); Vtot=m.d.var(ddof=0)
    Vbet=((m.groupby('document_chembl_id').d.transform('mean')-mu)**2).mean()
    c=m.groupby('molecule_chembl_id').filter(lambda x: x.document_chembl_id.nunique()>=2)
    out=dict(N=len(m),docs=m.document_chembl_id.nunique(),Vtot=Vtot,Vbet=Vbet,share=Vbet/Vtot,nrep=c.molecule_chembl_id.nunique())
    if out['nrep']<15: return out
    Vx=vcross(c); mols=c.molecule_chembl_id.unique(); idx={x:c[c.molecule_chembl_id==x] for x in mols}
    boots=[vcross(pd.concat([idx[x].assign(molecule_chembl_id=f'{x}_{i}') for i,x in enumerate(rng.choice(mols,len(mols)))])) for _ in range(nboot)]
    lo,hi=np.percentile(boots,[2.5,97.5])
    out.update(Vx=Vx,Vx_lo=lo,Vx_hi=hi,ceiling=1-Vx/Vtot,ceil_lo=1-hi/Vtot,ceil_hi=1-lo/Vtot,
        prot_share_of_between=min(Vx/Vbet,1.0) if Vbet>0 else np.nan,series_resid=max(Vbet-Vx,0)/Vtot)
    return out
if __name__=='__main__':
    os.makedirs('results',exist_ok=True)
    cache={}; rows=[]; cur=[]
    def get(g):
        t=T[g][0][0]
        if t not in cache:
            cache[t]=load(t)
            if not cache[t].empty:
                c=curate(cache[t]); c['target_chembl_id']=t; c['gene']=g; cur.append(c)
        return cache[t]
    for l in open('data/pairs.txt'):
        A,B=l.strip().split('/'); a,b=get(A),get(B)
        if a.empty or b.empty: continue
        s=pair_stats(a,b)
        if s: rows.append(dict(pair=f'{A}/{B}',**s)); print(rows[-1]['pair'],flush=True)
    pd.DataFrame(rows).to_csv('results/varcomp_by_pair.csv',index=False)
    pd.concat(cur).drop_duplicates(['molecule_chembl_id','document_chembl_id','target_chembl_id']).round(4)\
      .to_csv('data/curated_pX_by_compound_document_target.csv.gz',index=False)
