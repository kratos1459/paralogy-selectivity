#!/usr/bin/env python3
"""Stage 29: figures for the extension analyses (600 dpi PNG and LZW TIFF).

Figure E1  matched molecular pairs: offset removal, cross-publication reproducibility by effect size.
Figure E2  physicochemical drivers and pocket substitution chemistry.
Figure E3  representation benchmark (written only if stage 25 has finished).

usage: python 29_figures_ext.py
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common_ext import RES, FIG

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.linewidth": 1.0})
FAMC = {"kinase": "#1f4e9e", "epigenetic": "#7d3c98", "enzyme": "#1e8449", "gpcr": "#e67e22", "channel": "#b03a2e"}


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=600)
    fig.savefig(FIG / f"{name}.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw", "dpi": (600, 600)})
    plt.close(fig)
    print("wrote", name)


def lab(axes):
    for a, l in zip(np.ravel(axes), "abcdef"):
        a.text(-0.16, 1.05, f"({l})", transform=a.transAxes, fontsize=12, fontweight="bold")


# ---------------------------------------------------------------- Figure E1: MMPs
P = pd.read_csv(RES / "mmp_per_pair.csv")
ST = pd.read_csv(RES / "mmp_effect_size_strata.csv")
X = pd.read_csv(RES / "mmp_cross_document_effects.csv")
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
ax[0].scatter(P.between_share_selectivity, P.between_share_transformation,
              c=[FAMC.get(f, "#555") for f in P.family], s=26, edgecolor="k", lw=.4)
ax[0].plot([0, 1], [0, 1], "k--", lw=1)
ax[0].set_xlabel("Between-publication share,\nselectivity $\\Delta$")
ax[0].set_ylabel("Between-publication share,\nMMP transformation effect")
ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1)
ax[0].set_title("Matched pairs remove the offset", fontsize=9)
sub = X.sample(min(6000, len(X)), random_state=0)
ax[1].scatter(sub.effect_doc1, sub.effect_doc2, s=3, alpha=.25, color="#1f4e9e", edgecolor="none")
lim = 3
ax[1].plot([-lim, lim], [-lim, lim], "k--", lw=1); ax[1].axhline(0, color="k", lw=.6); ax[1].axvline(0, color="k", lw=.6)
ax[1].set_xlim(-lim, lim); ax[1].set_ylim(-lim, lim)
ax[1].set_xlabel("Transformation effect,\nfirst publication (log units)")
ax[1].set_ylabel("Transformation effect,\nsecond publication")
ax[1].set_title("Replication of substituent effects", fontsize=9)
b = ax[2].bar(range(len(ST)), ST.sign_concordance, color="#2e86c1", edgecolor="k", lw=.5)
ax[2].axhline(0.5, color="k", ls="--", lw=1)
ax2 = ax[2].twinx()
ax2.plot(range(len(ST)), ST.shrinkage, "o-", color="#b03a2e", ms=4, lw=1.4)
ax2.set_ylabel("Median |effect| ratio\n(second / first publication)", color="#b03a2e", fontsize=8)
ax2.tick_params(axis="y", colors="#b03a2e")
ax[2].set_xticks(range(len(ST))); ax[2].set_xticklabels(ST.bin, fontsize=8)
ax[2].set_xlabel("|Effect| in first publication (log units)")
ax[2].set_ylabel("Sign concordance"); ax[2].set_ylim(0.4, 0.8)
ax[2].set_title("Only large effects replicate", fontsize=9)
lab(ax)
save(fig, "FigureE1_matched_molecular_pairs")

# ---------------------------------------------------------------- Figure E2: physchem + pocket chemistry
S = pd.read_csv(RES / "physchem_summary.csv")
B = pd.read_csv(RES / "physchem_best_property.csv")
A = pd.read_csv(RES / "pocket_chemistry_associations.csv")
M = pd.read_csv(RES / "pocket_chemistry_merged.csv")
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4))
S = S.sort_values("median_abs_rho_within")
y = np.arange(len(S)); w = 0.4
ax[0].barh(y - w / 2, S.median_abs_rho_raw, w, color="#aec3df", edgecolor="k", lw=.4, label="raw")
ax[0].barh(y + w / 2, S.median_abs_rho_within, w, color="#1f4e9e", edgecolor="k", lw=.4, label="within publication")
ax[0].set_yticks(y); ax[0].set_yticklabels(S["prop"], fontsize=7.5)
ax[0].set_xlabel("Median |Spearman $\\rho$| with selectivity")
ax[0].legend(fontsize=7); ax[0].set_title("Bulk properties explain little", fontsize=9)
ax[1].hist(B.r2_best_property, bins=18, color="#1e8449", edgecolor="k", lw=.4)
ax[1].axvline(B.r2_best_property.median(), color="k", ls="--", lw=1.2)
ax[1].set_xlabel("$R^2$ of best single property\n(within-series selectivity)")
ax[1].set_ylabel("Target pairs"); ax[1].set_title("Within-series variance explained", fontsize=9)
sel = A[A.target == "sd_delta_samedoc"].sort_values("rho")
yy = np.arange(len(sel))
ax[2].barh(yy, sel.rho, color=["#b03a2e" if p >= 0.05 else "#1f4e9e" for p in sel.p_perm_within_family],
           edgecolor="k", lw=.4)
ax[2].errorbar(sel.rho, yy, xerr=[sel.rho - sel.ci_lo, sel.ci_hi - sel.rho], fmt="none", ecolor="k", lw=.8, capsize=2)
ax[2].axvline(0, color="k", lw=.8)
ax[2].set_yticks(yy); ax[2].set_yticklabels(sel.feature, fontsize=7.5)
ax[2].set_xlabel("Spearman $\\rho$ with selectivity window")
ax[2].set_title("Pocket substitution chemistry", fontsize=9)
lab(ax)
save(fig, "FigureE2_chemistry_determinants")

# ---------------------------------------------------------------- Figure E3: representations
f = RES / "representation_benchmark_per_pair.csv"
if f.exists():
    D = pd.read_csv(f)
    order = D.groupby("representation").r2.median().sort_values(ascending=False).index.tolist()
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.5))
    data = [D[D.representation == r].r2.values for r in order]
    bp = ax[0].boxplot(data, patch_artist=True, showfliers=False, widths=.6)
    for p in bp["boxes"]: p.set(facecolor="#aec3df", edgecolor="k", lw=.6)
    for m in bp["medians"]: m.set(color="k", lw=1.4)
    ax[0].axhline(0, color="k", ls="--", lw=1)
    ax[0].set_xticks(range(1, len(order) + 1)); ax[0].set_xticklabels(order, rotation=35, ha="right", fontsize=7)
    ax[0].set_ylabel("Leave-publication-out $R^2$")
    ax[0].set_title("Representation does not change the outcome", fontsize=9)
    dw = [D[D.representation == r].r2_within.values for r in order]
    bp2 = ax[1].boxplot(dw, patch_artist=True, showfliers=False, widths=.6)
    for p in bp2["boxes"]: p.set(facecolor="#cdebc5", edgecolor="k", lw=.6)
    for m in bp2["medians"]: m.set(color="k", lw=1.4)
    ax[1].axhline(0, color="k", ls="--", lw=1)
    ax[1].set_xticks(range(1, len(order) + 1)); ax[1].set_xticklabels(order, rotation=35, ha="right", fontsize=7)
    ax[1].set_ylabel("Within-series $R^2$")
    ax[1].set_title("Within-series ranking stays near zero", fontsize=9)
    lab(ax)
    save(fig, "FigureE3_representations")
else:
    print("stage 25 results not ready; skipped Figure E3")
