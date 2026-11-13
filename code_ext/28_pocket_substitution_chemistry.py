#!/usr/bin/env python3
"""Stage 28: what chemistry distinguishes the two binding pockets, and does it explain the
selectivity window?

The manuscript relates selectivity spread to a scalar pocket divergence (fraction of non-identical
pocket positions). That says how many residues differ, not how they differ. Here every differing
pocket position is classified chemically - change in side-chain volume, Kyte-Doolittle
hydropathy, formal charge, aromaticity and hydrogen-bonding capacity - and the resulting
chemistry-aware descriptors are tested against the empirical selectivity window (SD of
same-publication delta) and against leave-publication-out R2.

Inputs are the cached pocket residue positions (data/pdb/<ACC>.pocket.json) and UniProt sequences
(data/pdb/<ACC>.fa) produced by code/20_pocket_divergence.py; no network access is needed.

usage: python 28_pocket_substitution_chemistry.py
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from collections import Counter
from Bio import Align
from Bio.Align import substitution_matrices
from scipy.stats import spearmanr
from common_ext import ROOT, RES

PDB = ROOT / "data" / "pdb"
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2, I=166.7,
           L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8,
          K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
CHG = dict(D=-1, E=-1, K=1, R=1, H=0.1)
AROM = set("FWYH")
HBD = set("STYNQKRWH")
HBA = set("STYNQDEH")

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score, aligner.extend_gap_score, aligner.mode = -10, -0.5, "global"


def seq(acc):
    p = PDB / f"{acc}.fa"
    return "".join(open(p).read().split("\n")[1:]) if p.exists() else None


def pocket_positions(acc):
    p = PDB / f"{acc}.pocket.json"
    if not p.exists(): return None
    ent = json.load(open(p))
    if not ent: return None
    if len(ent) >= 3:
        c = Counter(r for e in ent for r in set(e["res"]))
        return {r for r, n in c.items() if n >= 2}
    return {r for e in ent for r in e["res"]}


def align_map(sA, sB):
    aln = aligner.align(sA, sB)[0]
    a, b = aln.aligned
    mapping = {}
    for (s1, e1), (s2, e2) in zip(a, b):
        for k in range(e1 - s1):
            mapping[s1 + k + 1] = s2 + k + 1          # 1-based UniProt positions
    return mapping


pairs = pd.read_csv(ROOT / "data" / "pairs.csv")
tbl = pd.read_csv(ROOT / "results" / "pocket_pair_table.csv")
rows, subs_rows = [], []
for pr in pairs.itertuples():
    sA, sB = seq(pr.accA), seq(pr.accB)
    pA, pB = pocket_positions(pr.accA), pocket_positions(pr.accB)
    if not (sA and sB and pA and pB): continue
    m = align_map(sA, sB)
    inv = {v: k for k, v in m.items()}
    cols = set()
    for r in pA:
        if r in m: cols.add((r, m[r]))
    for r in pB:
        if r in inv: cols.add((inv[r], r))
    if len(cols) < 5: continue
    subs = []
    for ra, rb in sorted(cols):
        if ra > len(sA) or rb > len(sB): continue
        a, b = sA[ra - 1], sB[rb - 1]
        if a not in VOL or b not in VOL: continue
        subs.append(dict(posA=ra, posB=rb, resA=a, resB=b, same=a == b,
                         dvol=VOL[b] - VOL[a], dkd=KD[b] - KD[a],
                         dchg=CHG.get(b, 0) - CHG.get(a, 0),
                         darom=int(b in AROM) - int(a in AROM),
                         dhbd=int(b in HBD) - int(a in HBD),
                         dhba=int(b in HBA) - int(a in HBA)))
    if not subs: continue
    S = pd.DataFrame(subs)
    diff = S[~S.same]
    name = f"{pr.A}/{pr.B}"
    rows.append(dict(pair=name, family=pr.family, n_pocket_positions=len(S), n_substitutions=len(diff),
                     pocket_divergence=len(diff) / len(S),
                     sum_abs_dvol=float(diff.dvol.abs().sum()), max_abs_dvol=float(diff.dvol.abs().max()) if len(diff) else 0.0,
                     mean_abs_dvol=float(diff.dvol.abs().mean()) if len(diff) else 0.0,
                     sum_abs_dkd=float(diff.dkd.abs().sum()), mean_abs_dkd=float(diff.dkd.abs().mean()) if len(diff) else 0.0,
                     n_charge_changes=int((diff.dchg.abs() > 0.5).sum()),
                     n_aromatic_changes=int((diff.darom != 0).sum()),
                     n_hbond_changes=int(((diff.dhbd != 0) | (diff.dhba != 0)).sum())))
    for r in diff.itertuples():
        subs_rows.append(dict(pair=name, family=pr.family, posA=r.posA, resA=r.resA, posB=r.posB, resB=r.resB,
                              dvol=r.dvol, dkd=r.dkd, dchg=r.dchg, darom=r.darom))
    print(f"{name}: {len(diff)}/{len(S)} pocket positions differ", flush=True)

D = pd.DataFrame(rows)
P = pd.DataFrame(subs_rows)
D.to_csv(RES / "pocket_chemistry_per_pair.csv", index=False)
P.to_csv(RES / "pocket_substitutions.csv", index=False)

key = tbl.assign(pair=tbl.A + "/" + tbl.B)[["pair", "sd_delta_samedoc", "r2_direct", "r2_within_direct", "pocket_div", "global_div", "frac_abs_delta_ge1"]]
M = D.merge(key, on="pair", how="inner")
M.to_csv(RES / "pocket_chemistry_merged.csv", index=False)

feats = ["pocket_divergence", "n_substitutions", "sum_abs_dvol", "max_abs_dvol", "mean_abs_dvol",
         "sum_abs_dkd", "mean_abs_dkd", "n_charge_changes", "n_aromatic_changes", "n_hbond_changes"]
out = []
rng = np.random.default_rng(0)
for f in feats:
    for tgt in ["sd_delta_samedoc", "frac_abs_delta_ge1", "r2_direct", "r2_within_direct"]:
        s = M[[f, tgt]].dropna()
        if len(s) < 10: continue
        rho = spearmanr(s[f], s[tgt])[0]
        bs = [spearmanr(s.sample(len(s), replace=True, random_state=int(rng.integers(1e9)))[f],
                        s.sample(len(s), replace=True, random_state=0)[tgt])[0] for _ in range(0)]
        boot = []
        for _ in range(2000):
            idx = rng.integers(0, len(s), len(s))
            boot.append(spearmanr(s[f].values[idx], s[tgt].values[idx])[0])
        # within-family permutation
        perm = []
        fam = M.loc[s.index, "family"].values
        for _ in range(5000):
            v = s[f].values.copy()
            for g in np.unique(fam):
                i = np.where(fam == g)[0]
                v[i] = rng.permutation(v[i])
            perm.append(spearmanr(v, s[tgt])[0])
        out.append(dict(feature=f, target=tgt, n=len(s), rho=float(rho),
                        ci_lo=float(np.nanpercentile(boot, 2.5)), ci_hi=float(np.nanpercentile(boot, 97.5)),
                        p_perm_within_family=float((np.sum(np.abs(perm) >= abs(rho)) + 1) / 5001)))
A = pd.DataFrame(out)
A.to_csv(RES / "pocket_chemistry_associations.csv", index=False)
print(A[A.target == "sd_delta_samedoc"].round(3).to_string(index=False))

top = P.assign(absv=P.dvol.abs()).sort_values("absv", ascending=False).head(20)
json.dump(dict(n_pairs=int(len(D)), median_substitutions=float(D.n_substitutions.median()),
               best_for_window=A[A.target == "sd_delta_samedoc"].sort_values("rho", ascending=False).head(3).to_dict("records"),
               best_for_lpo=A[A.target == "r2_direct"].sort_values("rho", ascending=False).head(3).to_dict("records"),
               largest_volume_substitutions=top[["pair", "resA", "posA", "resB", "posB", "dvol"]].to_dict("records")),
          open(RES / "pocket_chemistry_summary.json", "w"), indent=1, default=str)
print("wrote pocket chemistry results")
