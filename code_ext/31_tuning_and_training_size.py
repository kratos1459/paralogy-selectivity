#!/usr/bin/env python3
"""Stage 31: are the training-set caps and default hyperparameters responsible for the
leave-publication-out failure?

Limitation 5 of the manuscript concedes that training sets were capped at 2,000-6,000 compounds and
that hyperparameters were not tuned exhaustively. Two experiments answer it:

(a) training-size ablation - the same random forest at caps of 500, 1,000, 2,000, 4,000 and the full
    training set, under identical publication splits;
(b) nested tuning - random-forest and gradient-boosting hyperparameters selected inside each
    training set by publication-grouped 3-fold cross-validation, then applied once to the held-out
    publications, against the manuscript's defaults.

usage: python 31_tuning_and_training_size.py [start] [end] [tag]
"""
import sys, json, itertools
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from common_ext import load_all, matched, lpo_frame, lpo_splits, metrics, featurise, RES

try:
    import xgboost as xgb
    HAVE_XGB = True
except Exception:
    HAVE_XGB = False

NJOBS, NSPLIT, INNER_CAP = 4, 3, 1500
CAPS = [500, 1000, 2000, 4000, None]
RF_GRID = [dict(max_features="sqrt", min_samples_leaf=1), dict(max_features="sqrt", min_samples_leaf=5),
           dict(max_features=0.3, min_samples_leaf=1), dict(max_features=1.0, min_samples_leaf=5)]
XGB_GRID = [dict(max_depth=3, learning_rate=0.05), dict(max_depth=6, learning_rate=0.05),
            dict(max_depth=6, learning_rate=0.2)]


def rf(seed, **kw):
    return RandomForestRegressor(n_estimators=kw.pop("n_estimators", 150), n_jobs=NJOBS, random_state=seed, **kw)


def inner_docs_cv(X, y, docs, make, seed):
    """Publication-grouped 3-fold CV inside the training set; returns mean R2."""
    ud = np.array(sorted(set(docs)))
    rng = np.random.default_rng(seed); rng.shuffle(ud)
    folds = np.array_split(ud, 3)
    sc = []
    for f in folds:
        te = np.isin(docs, f)
        if te.sum() < 10 or (~te).sum() < 30: continue
        idx = np.where(~te)[0]
        if len(idx) > INNER_CAP: idx = rng.choice(idx, INNER_CAP, replace=False)
        m = make(seed).fit(X[idx], y[idx])
        sc.append(r2_score(y[te], m.predict(X[te])))
    return float(np.mean(sc)) if sc else -9


doc, smi, pairs = load_all()
start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(pairs)
tag = sys.argv[3] if len(sys.argv) > 3 else ""
size_rows, tune_rows = [], []
for i, pr in enumerate(pairs.itertuples()):
    if i < start or i >= end: continue
    mg = lpo_frame(matched(doc, pr.tidA, pr.tidB))
    if mg is None: continue
    name = f"{pr.A}/{pr.B}"
    X, keep = featurise(mg.molecule_chembl_id.values, smi, "morgan2_2048")
    if keep.sum() < 60: continue
    sub = mg[keep].reset_index(drop=True)
    y = sub.d.values; dd = sub.document_chembl_id.values

    # ---- (a) training-size ablation
    for cap in CAPS:
        got = []
        for s, tr, te in lpo_splits(sub, NSPLIT):
            idx = np.where(tr)[0]
            if cap and len(idx) > cap:
                idx = np.random.default_rng(s).choice(idx, cap, replace=False)
            m = rf(s, max_features="sqrt", min_samples_leaf=1, n_estimators=100).fit(X[idx], y[idx])
            got.append(dict(n_train=len(idx), **metrics(y[te], m.predict(X[te]), dd[te])))
        if got:
            size_rows.append(dict(pair=name, family=pr.family, cap=("all" if cap is None else cap),
                                  n_train=float(np.mean([g["n_train"] for g in got])),
                                  **{k: float(np.mean([g[k] for g in got])) for k in ("r2", "r2_within", "rho_within", "k")}))

    # ---- (b) nested tuning
    got_def, got_rf, got_xgb, chosen = [], [], [], []
    for s, tr, te in lpo_splits(sub, NSPLIT):
        idx = np.where(tr)[0]
        if len(idx) > 4000: idx = np.random.default_rng(s).choice(idx, 4000, replace=False)
        Xtr, ytr, dtr = X[idx], y[idx], dd[idx]
        got_def.append(metrics(y[te], rf(s, max_features="sqrt", min_samples_leaf=1, n_estimators=100).fit(Xtr, ytr).predict(X[te]), dd[te]))
        best, bestsc = None, -9
        for g in RF_GRID:
            sc = inner_docs_cv(Xtr, ytr, dtr, lambda sd, g=g: rf(sd, **dict(g)), s)
            if sc > bestsc: best, bestsc = g, sc
        chosen.append(str(best))
        got_rf.append(metrics(y[te], rf(s, **dict(best)).fit(Xtr, ytr).predict(X[te]), dd[te]))
        if HAVE_XGB:
            bx, bxs = None, -9
            for g in XGB_GRID:
                mk = lambda sd, g=g: xgb.XGBRegressor(n_estimators=200, subsample=0.8, colsample_bytree=0.5,
                                                      n_jobs=NJOBS, random_state=sd, tree_method="hist", **g)
                sc = inner_docs_cv(Xtr, ytr, dtr, mk, s)
                if sc > bxs: bx, bxs = g, sc
            mdl = xgb.XGBRegressor(n_estimators=200, subsample=0.8, colsample_bytree=0.5, n_jobs=NJOBS,
                                   random_state=s, tree_method="hist", **bx).fit(Xtr, ytr)
            got_xgb.append(metrics(y[te], mdl.predict(X[te]), dd[te]))
    if got_def:
        row = dict(pair=name, family=pr.family, n_rows=len(sub), n_docs=sub.document_chembl_id.nunique(),
                   chosen_rf=";".join(sorted(set(chosen))))
        for lbl, g in (("default", got_def), ("tuned_rf", got_rf), ("tuned_xgb", got_xgb)):
            if g:
                for k in ("r2", "r2_within", "rho_within"):
                    row[f"{lbl}_{k}"] = float(np.mean([x[k] for x in g]))
        tune_rows.append(row)
    print(f"{i+1} {name} done", flush=True)

pd.DataFrame(size_rows).to_csv(RES / f"training_size_ablation{tag}.csv", index=False)
pd.DataFrame(tune_rows).to_csv(RES / f"nested_tuning{tag}.csv", index=False)
print("wrote", tag or "(full)")
