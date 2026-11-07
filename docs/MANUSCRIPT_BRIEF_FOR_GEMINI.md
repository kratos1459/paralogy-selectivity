# PARALOG — manuscript assembly brief (target: J. Chem. Inf. Model., Full-Length Article)

**Role for the assembler:** build one coherent manuscript from the blocks below. Do not invent numbers, citations or analyses; every number here is from the released pipeline (https://github.com/kratos1459/paralogy-selectivity, `results/`). Keep ACS style. Mark any gap as `[TODO]` instead of filling it.

## 0. Metadata
- **Title:** Evaluation Protocol Limits Learnability of Paralog Selectivity Across ChEMBL Target Pairs
- **Authors:** Tawsif Al Arian (corresponding; tawsifalarian43@gmail.com), Most. Afrin Akter, Md. Shoriful Alam Sakib, Ehsanul Islam (co-senior contact; ehsanul.islam@juniv.edu). Affiliations: Dept. of Pharmacy (Al Arian, Akter, Islam) and Dept. of Biotechnology and Genetic Engineering (Sakib), Jahangirnagar University, Savar, Dhaka-1342, Bangladesh.
- **Keywords:** paralog selectivity; isoform selectivity; ChEMBL; QSAR; scaffold split; leave-publication-out; noise ceiling; few-shot calibration.
- **One-sentence story (keep it the spine):** Standard splits overstate selectivity-model accuracy because most selectivity variance sits *between publications*; we quantify this (variance decomposition, noise ceiling), show it on hand-picked and random pairs, and show that one or two measured anchor compounds with empirical-Bayes shrinkage recover much of the lost accuracy.

## 1. Data and definitions (Methods)
- ChEMBL 37; human single-protein targets; IC50/Ki/Kd in nM; relation "="; no data_validity_comment; potential_duplicate = 0. pX = 9 − log10(nM), 2 ≤ pX ≤ 12.5. Median-aggregated per (compound, document, target).
- Selectivity Δ = pX(A) − pX(B) for a compound measured on both targets in the **same document**.
- **Curated set:** 51 hand-selected paralog/isoform pairs (93 targets; kinases 15, epigenetic 6, other enzymes 18, GPCR 10, ion channels 2). 47 pairs have enough data for model evaluation.
- **Random set (new):** 60 pairs drawn at random from ChEMBL protein-classification leaf classes (class level ≥3), human single-protein targets with ≥500 exact records, ≤3 pairs per class, original 51 pairs excluded, seed 1 (`code/17_sample_random_pairs.py`). 45 pass filters (≥60 matched compounds); 35 have enough multi-compound documents (≥5 matched compounds, ≥8 such documents, ≥60 compounds) for leave-publication-out (LPO).
- Evaluation regimes: Bemis–Murcko scaffold split (25 % of scaffold groups, 5 seeds); LPO (25 % of documents held out, test compounds removed from training); chronological split (median-year cut, 2017).
- Model: random forest, Morgan count fingerprints (2048 bits, radius 2). Comparators: HistGradientBoosting, MLP, XGBoost 3.4.1, series-centred RF.
- Variance decomposition: Vtot, between-publication share, within-compound cross-document variance Vx, **noise ceiling = 1 − Vx/Vtot** (500× bootstrap).
- **Anchor calibration:** m measured compounds in the new publication shift predictions by an offset. Plain: mean residual. Empirical Bayes: μ + w(r − μ), w = τ²/(τ² + σ²/m), τ² and σ² from publication-grouped nested CV (out-of-bag residuals underestimate τ²: reported negative result). Theory: MSE_plain = σ²(1+1/m), MSE_EB = σ²(1+w/m) (`docs/THEORY.md`).

## 2. Results blocks (numbers are medians over pairs unless stated)
**R1 Scaffold vs LPO (curated).** Scaffold R² 0.54; LPO 0.26; within-series 0.04; scaffold leakage 92 %; between-publication variance share 65 % (IQR 54–76 %).
**R2 Noise ceiling (curated, 30 pairs with ≥20 repeat compounds).** Median ceiling 0.63 (IQR 0.37–0.71). LPO R² = 55 % of ceiling and never exceeds it; scaffold R² exceeds the ceiling in 6/26 pairs (leakage signature). Ceiling predicts LPO R² (Spearman 0.73, p = 5×10⁻⁶). ≈57 % of between-publication variance is protocol/noise; series component ≈27 % of total variance. Chronological split median R² 0.07 (45 pairs), within-series 0.02.
**R3 Anchor calibration (curated, 47 pairs, 4 fresh splits).** Plain vs EB pooled R²: m=1 0.18→0.34 (EB better in 46/47, p=2×10⁻¹²); m=2 0.38→0.46 (41/47); m=3 0.44→0.47 (41/47); m=5 0.49→0.50 (36/47). One EB anchor ≈ two plain anchors. Theory vs observed MSE: Pearson 0.95–0.99 (gain slightly overestimated, 0.137 vs 0.117 at m=1).
**R4 Learners.** RF 0.27, HGB 0.17, MLP 0.17 pooled R² (cap 2500); XGBoost 0.14 vs RF 0.20 (RF better 36/47, p=2×10⁻⁵). Stronger learners do not rescue cross-publication prediction.
**R5 Assay descriptors / assay type.** Recorded assay descriptors explain ≈20 % of cross-document variance (pooled 0.57→0.46; median 34 % per pair). Between-publication share 0.56 for Ki/Kd (37 pairs) vs 0.66 for IC50 (45 pairs).
**R6 Single-document hold-out (honest limitation).** 47 tests/34 documents: base R² median −0.11; EB beats plain in 42–44/47 but R² stays ≈ ≤0.
**R7 Random-pair replication (NEW; `results/random_pairs/`, Figure 7).** 45 pairs pass filters; 35 for LPO. Scaffold R² 0.55 (curated 0.54); **LPO R² −0.10** (curated 0.26); scaffold leakage 0.91; between-publication share 0.70 (curated 0.68); ceiling 0.65. EB vs plain (LPO R², m=1/2/3/5): EB 0.31/0.37/0.41/0.37 vs plain 0.12/0.32/0.38/0.40 → clear gain at m=1–2, none at m=5. Interpretation: the diagnosis (leakage, between-publication variance) is not an artefact of pair selection; absolute LPO accuracy is lower on random pairs because they have fewer/smaller publications; no unanchored model predicts across publications there.
*Caveat to state:* curated vs random LPO is not a like-for-like difficulty comparison (different document counts); report n for each.

## 3. Figures and tables (files in the repo)
- Fig 1–5: existing (see `Manuscript_JCIM_ACS.docx`).
- **Fig 6** `results/figure_ceiling_and_calibration.png` — (a) R² vs noise ceiling per pair, scaffold vs LPO; (b) plain vs EB pooled R² against anchors m.
- **Fig 7** `results/figure_random_pair_replication.png` — (a) scaffold vs LPO R² for curated vs random pairs; (b) between-publication share histograms; (c) anchor curves, plain (dashed) vs EB (solid). Palette: blue = curated, orange = random.
- **Table 5** plain vs EB at m = 0,1,2,3,5 (values in R3). Suggested new **Table 6**: curated vs random summary (rows: n pairs, scaffold R², LPO R², leakage, between-publication share, ceiling, EB gain at m=1).
- TOC graphic (3.25 × 1.75 in): `[TODO]` — suggest scaffold-split vs LPO bar pair plus an "anchor" arrow.

## 4. Suggested structure and framing for a JCIM paper
1. Abstract (≤250 words; lead with the audit and the ceiling; state random-pair replication and the m=1 EB result; state the limitation that unanchored cross-publication accuracy is low).
2. Introduction — gap: selectivity benchmarks mix publications; no quantification of how much of Δ variance is protocol/series.
3. Methods — §1 above, one pipeline throughout (**use the released re-implementation's numbers everywhere**; the original manuscript numbers differ slightly: LPO R² 0.21 vs 0.26).
4. Results — R1 → R2 → R7 → R3 → R4/R5 → R6 (limitation last, not hidden).
5. Discussion — implications: report LPO + ceiling; use anchors; ceiling is estimated from repeat (mostly reference) compounds and may be biased; no prospective/external validation.
6. Conclusions; Data and Software Availability (GitHub, MIT; Zenodo DOI https://doi.org/10.5281/zenodo.24035853); Supporting Information list; Author contributions `[TODO]`; Funding/competing interests `[TODO]`.
- Do **not** claim state-of-the-art prediction. Claim: a quantified evaluation audit, a noise-ceiling diagnostic, a calibrated few-shot fix.

## 5. Known limitations the text must keep
Single database (ChEMBL); median-aggregation of repeated measurements; noise ceiling needs repeat compounds; random pairs: 35/60 usable for LPO; no pretrained/graph models; no external or prospective validation; EB gain vanishes by m=5 on random pairs.
