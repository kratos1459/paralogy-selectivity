#!/usr/bin/env python3
"""Stage 26: do substituent changes carry reproducible selectivity effects across publications?

The whole-molecule selectivity endpoint is dominated by a publication offset, which is why
leave-publication-out models fail. A matched molecular pair (MMP) measures the effect of a single
substituent change *inside* one publication, so the offset cancels by construction. This asks the
chemistry question directly: is the structure-selectivity relationship of a given transformation
(for example para-Cl -> para-OMe) reproducible between independent laboratories, even though
whole-molecule models are not?

For every target pair:
 - molecules are cut at each acyclic single bond between heavy atoms, giving (core, R-group) pairs;
 - within one publication, two molecules sharing a core form an MMP with effect dd = delta_i - delta_j
   for the transformation T = (R_j -> R_i);
 - transformation effects are aggregated per publication, and the between-publication variance share
   of those effects is compared with that of the raw selectivity values;
 - a transfer model predicts the effect of each MMP in a held-out publication from the median effect
   of the same transformation in training publications, and is scored against predicting zero.

usage: python 26_matched_molecular_pairs.py [n_pairs]
"""
import sys, json
from collections import defaultdict
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sklearn.metrics import r2_score
from scipy.stats import spearmanr, wilcoxon
from common_ext import load_all, matched, lpo_frame, RES, _mol

MAX_HEAVY, MIN_CORE_ATOMS, MAX_R_ATOMS = 70, 6, 20


def cuts(mol):
    """All single acyclic C-C/C-X bond cuts -> list of (core_smiles, r_smiles)."""
    out = []
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing(): continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1: continue
        try:
            frag = Chem.FragmentOnBonds(mol, [b.GetIdx()], addDummies=True, dummyLabels=[(1, 1)])
            pieces = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=True)
        except Exception:
            continue
        if len(pieces) != 2: continue
        n = [p.GetNumHeavyAtoms() - 1 for p in pieces]
        big, small = (0, 1) if n[0] >= n[1] else (1, 0)
        if n[big] < MIN_CORE_ATOMS or n[small] > MAX_R_ATOMS or n[small] < 1: continue
        try:
            out.append((Chem.MolToSmiles(pieces[big]), Chem.MolToSmiles(pieces[small])))
        except Exception:
            continue
    return out


doc, smi, pairs = load_all()
limit = int(sys.argv[1]) if len(sys.argv) > 1 else len(pairs)
rows, pair_rows, trans_rows, cross = [], [], [], []
for i, pr in enumerate(pairs.itertuples()):
    if i >= limit: break
    mg = lpo_frame(matched(doc, pr.tidA, pr.tidB))
    if mg is None: continue
    name = f"{pr.A}/{pr.B}"
    # fragment every molecule once
    frag = {}
    for mid in mg.molecule_chembl_id.unique():
        m = _mol(mid, smi)
        if m is None or m.GetNumHeavyAtoms() > MAX_HEAVY: continue
        frag[mid] = cuts(m)
    if len(frag) < 20: continue

    # MMPs inside each publication
    recs = []
    for d, sub in mg.groupby("document_chembl_id"):
        dmap = dict(zip(sub.molecule_chembl_id, sub.d))
        core_idx = defaultdict(list)
        for mid in sub.molecule_chembl_id.unique():
            for core, r in frag.get(mid, []):
                core_idx[core].append((mid, r))
        for core, members in core_idx.items():
            seen = {}
            for mid, r in members:
                seen.setdefault((mid, r), True)
            uniq = list({(mid, r) for mid, r in members})
            if len(uniq) < 2: continue
            for a in range(len(uniq)):
                for b in range(a + 1, len(uniq)):
                    (m1, r1), (m2, r2) = uniq[a], uniq[b]
                    if m1 == m2 or r1 == r2: continue
                    dd = dmap[m2] - dmap[m1]
                    t = (r1, r2) if r1 < r2 else (r2, r1)
                    sign = 1.0 if (r1, r2) == t else -1.0
                    recs.append((d, core, f"{t[0]}>>{t[1]}", sign * dd))
    if len(recs) < 50: continue
    R = pd.DataFrame(recs, columns=["doc", "core", "transformation", "dd"])
    # one effect per (transformation, document)
    TD = R.groupby(["transformation", "doc"]).dd.median().reset_index()
    counts = TD.groupby("transformation").doc.nunique()
    multi = counts[counts >= 2].index
    if len(multi) < 10: continue
    M = TD[TD.transformation.isin(multi)]

    # variance decomposition of transformation effects vs raw selectivity
    mu = M.dd.mean()
    v_tot_t = ((M.dd - mu) ** 2).sum()
    v_bet_t = ((M.groupby("doc").dd.transform("mean") - mu) ** 2).sum()
    mu2 = mg.d.mean()
    v_tot_y = ((mg.d - mu2) ** 2).sum()
    v_bet_y = ((mg.groupby("document_chembl_id").d.transform("mean") - mu2) ** 2).sum()

    # cross-publication reproducibility of a transformation: first vs second document
    g = M.groupby("transformation").dd.apply(list)
    first = np.array([v[0] for v in g if len(v) >= 2]); second = np.array([v[1] for v in g if len(v) >= 2])
    rho = spearmanr(first, second)[0] if len(first) >= 8 else np.nan
    conc = float((np.sign(first) == np.sign(second)).mean()) if len(first) else np.nan
    for f_, s_ in zip(first, second):
        cross.append(dict(pair=name, family=pr.family, effect_doc1=float(f_), effect_doc2=float(s_)))

    # transfer model: predict each held-out document's MMP effects from training documents
    udocs = sorted(M.doc.unique())
    preds, obs = [], []
    for held in udocs:
        tr = M[M.doc != held]; te = M[M.doc == held]
        lib = tr.groupby("transformation").dd.median()
        hit = te[te.transformation.isin(lib.index)]
        if hit.empty: continue
        preds.extend(lib.loc[hit.transformation].values); obs.extend(hit.dd.values)
    if len(obs) < 30: continue
    obs = np.array(obs); preds = np.array(preds)
    r2_transfer = r2_score(obs, preds)
    r2_zero = r2_score(obs, np.zeros_like(obs))
    sign_acc = float((np.sign(obs) == np.sign(preds))[obs != 0].mean())

    pair_rows.append(dict(pair=name, family=pr.family, n_mmp=len(R), n_transformations=int(counts.size),
                          n_transformations_multi_doc=int(len(multi)), n_docs=int(mg.document_chembl_id.nunique()),
                          between_share_transformation=float(v_bet_t / v_tot_t) if v_tot_t > 0 else np.nan,
                          between_share_selectivity=float(v_bet_y / v_tot_y) if v_tot_y > 0 else np.nan,
                          sd_transformation_effect=float(M.dd.std()), sd_selectivity=float(mg.d.std()),
                          cross_doc_rho=float(rho) if rho == rho else np.nan, cross_doc_sign_concordance=conc,
                          n_transfer_predictions=int(len(obs)), r2_transfer=float(r2_transfer),
                          r2_predict_zero=float(r2_zero), sign_accuracy=sign_acc))
    top = M.groupby("transformation").agg(n_docs=("doc", "nunique"), median_effect=("dd", "median"),
                                          sd_effect=("dd", "std")).sort_values("n_docs", ascending=False).head(15)
    for t, r in top.iterrows():
        trans_rows.append(dict(pair=name, transformation=t, n_docs=int(r.n_docs),
                               median_effect=float(r.median_effect), sd_effect=float(r.sd_effect)))
    print(f"{i+1} {name}: {len(R)} MMPs, {len(multi)} transformations in >=2 docs, "
          f"between-share {v_bet_t/v_tot_t:.2f} vs {v_bet_y/v_tot_y:.2f}, transfer R2 {r2_transfer:.3f}", flush=True)

X = pd.DataFrame(cross)
if len(X):
    X.to_csv(RES / "mmp_cross_document_effects.csv", index=False)
    X["abs1"] = X.effect_doc1.abs()
    bins = [0, 0.3, 0.5, 1.0, 2.0, 99]
    X["bin"] = pd.cut(X.abs1, bins, labels=["<0.3", "0.3-0.5", "0.5-1.0", "1.0-2.0", ">2.0"])
    ST = X.groupby("bin", observed=True).apply(lambda g: pd.Series(dict(
        n=len(g), sign_concordance=float((np.sign(g.effect_doc1) == np.sign(g.effect_doc2)).mean()),
        median_abs_effect_doc1=float(g.abs1.median()), median_abs_effect_doc2=float(g.effect_doc2.abs().median()),
        spearman=float(spearmanr(g.effect_doc1, g.effect_doc2)[0]) if len(g) >= 8 else np.nan,
        shrinkage=float(np.median(np.abs(g.effect_doc2)) / max(np.median(np.abs(g.effect_doc1)), 1e-9)))), include_groups=False)
    ST.to_csv(RES / "mmp_effect_size_strata.csv")
    print(ST.round(3).to_string())

P = pd.DataFrame(pair_rows)
P.to_csv(RES / "mmp_per_pair.csv", index=False)
pd.DataFrame(trans_rows).to_csv(RES / "mmp_top_transformations.csv", index=False)
if len(P):
    s = dict(n_pairs=int(len(P)),
             median_between_share_transformation=float(P.between_share_transformation.median()),
             median_between_share_selectivity=float(P.between_share_selectivity.median()),
             wilcoxon_p=float(wilcoxon(P.between_share_transformation, P.between_share_selectivity).pvalue) if len(P) >= 6 else None,
             lower_in_pairs=int((P.between_share_transformation < P.between_share_selectivity).sum()),
             median_cross_doc_rho=float(P.cross_doc_rho.median()),
             median_sign_concordance=float(P.cross_doc_sign_concordance.median()),
             median_r2_transfer=float(P.r2_transfer.median()),
             median_sign_accuracy=float(P.sign_accuracy.median()),
             pairs_with_positive_transfer=int((P.r2_transfer > 0).sum()),
             median_sd_transformation=float(P.sd_transformation_effect.median()),
             median_sd_selectivity=float(P.sd_selectivity.median()),
             strata=ST.round(3).reset_index().astype(str).to_dict("records") if len(X) else None,
             overall_cross_doc_spearman=float(spearmanr(X.effect_doc1, X.effect_doc2)[0]) if len(X) > 8 else None,
             overall_sign_concordance=float((np.sign(X.effect_doc1) == np.sign(X.effect_doc2)).mean()) if len(X) else None)
    json.dump(s, open(RES / "mmp_summary.json", "w"), indent=1)
    print(json.dumps(s, indent=1))
