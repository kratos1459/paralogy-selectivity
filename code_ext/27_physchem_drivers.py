#!/usr/bin/env python3
"""Stage 27: which physicochemical properties drive paralog selectivity, once the publication
offset is removed?

Selectivity is correlated with a property both because the property genuinely shifts the balance
between two binding sites and because whole publications differ. Centring both the selectivity
endpoint and the property within each publication removes the offset, leaving the within-series
chemistry. For every target pair and every property we report the within-publication Spearman
correlation with selectivity, a within-publication permutation p-value, and the variance of
within-series selectivity explained by the single best property.

usage: python 27_physchem_drivers.py [n_pairs]
"""
import sys, json
import numpy as np, pandas as pd
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from common_ext import load_all, matched, lpo_frame, centered, RES, _mol

PROPS = {
    "MolWt": Descriptors.MolWt, "cLogP": Crippen.MolLogP, "TPSA": rdMolDescriptors.CalcTPSA,
    "HBD": rdMolDescriptors.CalcNumHBD, "HBA": rdMolDescriptors.CalcNumHBA,
    "RotatableBonds": rdMolDescriptors.CalcNumRotatableBonds,
    "AromaticRings": rdMolDescriptors.CalcNumAromaticRings,
    "AliphaticRings": rdMolDescriptors.CalcNumAliphaticRings,
    "FractionCSP3": rdMolDescriptors.CalcFractionCSP3,
    "HeavyAtoms": lambda m: m.GetNumHeavyAtoms(),
    "FormalCharge": lambda m: sum(a.GetFormalCharge() for a in m.GetAtoms()),
    "MolarRefractivity": Crippen.MolMR,
}
NPERM = 500


def perm_p(yc, xc, docs, obs, rng):
    """Permute the property within each publication."""
    hits = 0
    idx = pd.Series(range(len(xc))).groupby(np.asarray(docs))
    groups = [np.array(g) for _, g in idx]
    for _ in range(NPERM):
        xp = xc.copy()
        for g in groups:
            xp[g] = rng.permutation(xc[g])
        if abs(spearmanr(yc, xp)[0]) >= abs(obs): hits += 1
    return (1 + hits) / (NPERM + 1)


doc, smi, pairs = load_all()
limit = int(sys.argv[1]) if len(sys.argv) > 1 else len(pairs)
rows, best_rows = [], []
for i, pr in enumerate(pairs.itertuples()):
    if i >= limit: break
    mg = lpo_frame(matched(doc, pr.tidA, pr.tidB))
    if mg is None: continue
    name = f"{pr.A}/{pr.B}"
    mols = {m: _mol(m, smi) for m in mg.molecule_chembl_id.unique()}
    keep = mg.molecule_chembl_id.map(lambda m: mols.get(m) is not None).values
    sub = mg[keep].reset_index(drop=True)
    if len(sub) < 60: continue
    P = {}
    for pname, fn in PROPS.items():
        try:
            P[pname] = sub.molecule_chembl_id.map(lambda m: float(fn(mols[m]))).values
        except Exception:
            continue
    docs = sub.document_chembl_id.values
    yc = centered(sub.d.values, docs)
    rng = np.random.default_rng(0)
    feats = []
    for pname, v in P.items():
        if np.nanstd(v) == 0: continue
        xc = centered(v, docs)
        if np.nanstd(xc) == 0: continue
        rho_w = spearmanr(yc, xc)[0]
        rho_raw = spearmanr(sub.d.values, v)[0]
        p = perm_p(yc, xc, docs, rho_w, rng)
        rows.append(dict(pair=name, family=pr.family, prop=pname, n=len(sub),
                         rho_within=float(rho_w), rho_raw=float(rho_raw), perm_p=float(p)))
        feats.append((pname, xc))
    if feats:
        X = np.column_stack([f[1] for f in feats])
        lin = LinearRegression().fit(X, yc)
        r2_all = r2_score(yc, lin.predict(X))
        bestp, bestr2 = None, -9
        for pname, xc in feats:
            r2 = r2_score(yc, LinearRegression().fit(xc.reshape(-1, 1), yc).predict(xc.reshape(-1, 1)))
            if r2 > bestr2: bestp, bestr2 = pname, r2
        best_rows.append(dict(pair=name, family=pr.family, n=len(sub), best_property=bestp,
                              r2_best_property=float(bestr2), r2_all_properties=float(r2_all),
                              sd_within_series=float(np.std(yc))))
    print(f"{i+1} {name}", flush=True)

D = pd.DataFrame(rows); D.to_csv(RES / "physchem_per_pair.csv", index=False)
B = pd.DataFrame(best_rows); B.to_csv(RES / "physchem_best_property.csv", index=False)
S = D.groupby("prop").agg(pairs=("pair", "nunique"),
                          median_abs_rho_within=("rho_within", lambda s: float(np.median(np.abs(s)))),
                          median_abs_rho_raw=("rho_raw", lambda s: float(np.median(np.abs(s)))),
                          pairs_perm_p05=("perm_p", lambda s: int((s < 0.05).sum())),
                          pairs_abs_rho_gt_02=("rho_within", lambda s: int((np.abs(s) > 0.2).sum()))
                          ).sort_values("median_abs_rho_within", ascending=False)
S.to_csv(RES / "physchem_summary.csv")
print(S.round(3).to_string())
print(B.groupby("best_property").size().sort_values(ascending=False).to_string())
json.dump(dict(summary=S.round(4).reset_index().to_dict("records"),
               median_r2_best_property=float(B.r2_best_property.median()) if len(B) else None,
               median_r2_all_properties=float(B.r2_all_properties.median()) if len(B) else None,
               best_property_counts=B.groupby("best_property").size().to_dict() if len(B) else {}),
          open(RES / "physchem_summary.json", "w"), indent=1)
