"""Relate model performance to the cross-publication noise ceiling and summarise reproduction (run from repo root)."""
import pandas as pd,json
from scipy.stats import spearmanr
m=pd.read_csv('results/models_by_pair.csv'); v=pd.read_csv('results/varcomp_by_pair.csv'); t=pd.read_csv('results/time_split_by_pair.csv')
x=m.merge(v[['pair','ceiling','nrep']],on='pair'); x=x[(x.nrep>=20)]
e=x[x.r2_lpo.notna()&(x.ceiling>0.05)].copy(); e['frac_of_ceiling']=e.r2_lpo/e.ceiling
e.to_csv('results/ceiling_vs_models.csv',index=False)
rho=spearmanr(x.ceiling,x.r2_lpo,nan_policy='omit')
s=dict(scaffold_R2_median=m.r2_scaffold.median(),scaffold_pairs=int(m.r2_scaffold.notna().sum()),lpo_R2_median=m.r2_lpo.median(),lpo_pairs=int(m.r2_lpo.notna().sum()),
 within_series_R2_median=m.r2_within.median(),between_pub_share_median=m.frac_between.median(),between_pub_share_IQR=m.frac_between.quantile([.25,.75]).tolist(),
 scaffold_leakage_median=m.leak.median(),**{f'anchor{k}_R2_median':m[f'r2_anchor{k}'].median() for k in (0,1,2,3,5)},
 time_split_pairs=len(t),time_split_R2_median=t.r2_time.median(),time_split_within_R2_median=t.r2_within_time.median(),
 ceiling_pairs=len(e),ceiling_median=e.ceiling.median(),lpo_over_ceiling_median=e.frac_of_ceiling.median(),
 pairs_scaffold_R2_above_ceiling=int((e.r2_scaffold>e.ceiling).sum()),pairs_lpo_R2_above_ceiling=int((e.r2_lpo>e.ceiling).sum()),
 spearman_ceiling_vs_lpo_R2=rho.statistic,spearman_p=rho.pvalue)
json.dump(s,open('results/model_summary.json','w'),indent=1,default=float); print(json.dumps(s,indent=1,default=float))
