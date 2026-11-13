#!/usr/bin/env python3
"""Stage 25: does the molecular representation change cross-publication generalisation?

The manuscript evaluates one representation (2,048-bit Morgan count fingerprints, radius 2). A
reviewer will ask whether a richer or different representation closes the scaffold-to-LPO gap.
Seven representations are run through the identical leave-publication-out protocol with the same
random forest: Morgan radius 2 and 3, atom pairs, topological torsions, MACCS keys, 208 RDKit
physicochemical descriptors, and Morgan concatenated with those descriptors.

usage: python 25_representations.py [n_pairs]
"""
import sys, json
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor
from common_ext import load_all, matched, lpo_frame, lpo_splits, metrics, featurise, REPS, RES

NSPLIT, CAP = 3, 4000
doc, smi, pairs = load_all()
start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
limit = int(sys.argv[2]) if len(sys.argv) > 2 else len(pairs)
tag = sys.argv[3] if len(sys.argv) > 3 else ""
rows = []
for i, pr in enumerate(pairs.itertuples()):
    if i < start or i >= limit: continue
    m = matched(doc, pr.tidA, pr.tidB)
    mg = lpo_frame(m)
    if mg is None: continue
    name = f"{pr.A}/{pr.B}"
    mids = mg.molecule_chembl_id.values
    for rep in REPS:
        X, keep = featurise(mids, smi, rep)
        if keep.sum() < 60: continue
        sub = mg[keep].reset_index(drop=True)
        y = sub.d.values; dd = sub.document_chembl_id.values
        got = []
        for s, tr, te in lpo_splits(sub, NSPLIT):
            idx = np.where(tr)[0]
            if len(idx) > CAP:
                idx = np.random.default_rng(s).choice(idx, CAP, replace=False)
            mod = RandomForestRegressor(100, max_features="sqrt", n_jobs=5, random_state=s).fit(X[idx], y[idx])
            got.append(metrics(y[te], mod.predict(X[te]), dd[te]))
        if not got: continue
        rows.append(dict(pair=name, family=pr.family, representation=rep, n_features=X.shape[1],
                         n_rows=len(sub), n_docs=sub.document_chembl_id.nunique(), n_splits=len(got),
                         **{k: float(np.mean([g[k] for g in got])) for k in got[0]}))
    print(f"{i+1}/{len(pairs)} {name}", flush=True)

D = pd.DataFrame(rows)
D.to_csv(RES / f"representation_benchmark_per_pair{tag}.csv", index=False)
if tag:
    raise SystemExit(0)
S = D.groupby("representation").agg(pairs=("pair", "nunique"), n_features=("n_features", "max"),
                                    lpo_r2=("r2", "median"), lpo_r2_iqr_lo=("r2", lambda s: s.quantile(.25)),
                                    lpo_r2_iqr_hi=("r2", lambda s: s.quantile(.75)),
                                    within_r2=("r2_within", "median"), within_rho=("rho_within", "median"),
                                    k=("k", "median"), r_centred=("r_centred", "median")).sort_values("lpo_r2", ascending=False)
S.to_csv(RES / "representation_benchmark_summary.csv")
print(S.round(3).to_string())

# paired comparison of every representation against the Morgan baseline
from scipy.stats import wilcoxon
base = D[D.representation == "morgan2_2048"].set_index("pair")
comp = []
for rep in REPS:
    if rep == "morgan2_2048": continue
    o = D[D.representation == rep].set_index("pair")
    j = base.join(o, lsuffix="_base", rsuffix="_rep", how="inner")
    if len(j) < 5: continue
    comp.append(dict(representation=rep, n_pairs=len(j),
                     median_delta_lpo_r2=float((j.r2_rep - j.r2_base).median()),
                     better_than_morgan=int((j.r2_rep > j.r2_base).sum()),
                     wilcoxon_p=float(wilcoxon(j.r2_rep, j.r2_base).pvalue),
                     median_delta_within_r2=float((j.r2_within_rep - j.r2_within_base).median())))
C = pd.DataFrame(comp); C.to_csv(RES / "representation_vs_morgan.csv", index=False)
print(C.round(4).to_string(index=False))
json.dump(dict(summary=S.round(4).reset_index().to_dict("records"), vs_morgan=C.round(4).to_dict("records")),
          open(RES / "representation_summary.json", "w"), indent=1)
