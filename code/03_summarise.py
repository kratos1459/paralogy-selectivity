"""Summarise results/varcomp_by_pair.csv -> results/summary.json (run from repo root)."""
import pandas as pd,json
d=pd.read_csv('results/varcomp_by_pair.csv'); e=d[d.nrep>=20].copy(); e['ratio']=e.Vx/e.Vtot
s=dict(pairs_analysed=len(d),between_publication_share_median=d.share.median(),
 between_publication_share_IQR=d.share.quantile([.25,.75]).tolist(),pairs_with_ge20_repeat_compounds=len(e),
 noise_ratio_Vx_over_Vtot_median=e.ratio.median(),noise_ratio_IQR=e.ratio.quantile([.25,.75]).tolist(),
 noise_ceiling_R2_median=e.ceiling.median(),noise_ceiling_R2_IQR=e.ceiling.quantile([.25,.75]).tolist(),
 pairs_with_nonpositive_ceiling=int((e.ceiling<=0).sum()),
 protocol_share_of_between_median=e.prot_share_of_between.median(),series_share_of_total_median=e.series_resid.median())
json.dump(s,open('results/summary.json','w'),indent=1,default=float); print(json.dumps(s,indent=1,default=float))
