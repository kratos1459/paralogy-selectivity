# PARALOG: variance components of paralog-selectivity data

Code and derived data for the variance-component analysis that accompanies the manuscript on
leave-publication-out (LPO) evaluation of paralog/isoform selectivity models (ChEMBL 37, 51 target pairs).

## What this repository reproduces

The original 01–16 script pipeline lived in a sandbox that no longer exists, so the model analyses were **re-implemented** from the manuscript's Methods and re-run on a fresh ChEMBL 37 pull. The re-implementation matches the manuscript closely (below). Not re-implemented: XGBoost and graph-network comparison (scikit-learn gradient boosting and MLP used instead), pairwise training, and the PDB pocket-divergence analysis.

| Quantity | Manuscript | This repository |
|---|---|---|
| Scaffold-split R², median (pairs) | 0.56 (49) | 0.54 (49) |
| Leave-publication-out R², median (pairs) | 0.21 (47) | 0.26 (47) |
| Within-series R², median | 0.05 | 0.04 |
| Between-publication variance share, median (IQR) | 65% (54–76%) | 65% (54–76%) |
| Scaffold-test compounds whose publication is in training | 92% | 92% |
| Anchor calibration, pooled R² at 0/1/2/3/5 anchors | 0.19/0.20/0.42/0.48/0.52 | 0.26/0.27/0.40/0.46/0.51 |

Remaining differences (largest: leave-publication-out R² 0.26 vs 0.21) come from unrecoverable details such as random seeds, hyperparameters and pair-selection thresholds; the qualitative picture is unchanged.

## New analyses added here

1. **Noise ceiling** (`02`, `06`). The same compound's selectivity varies across publications; this cross-publication replicate variance (protocol offset plus measurement error) bounds the R² any model can reach on a new publication: ceiling = 1 − Vx/Vtot. Median ceiling 0.63 (IQR 0.37 to 0.71; 30 pairs with at least 20 repeat compounds; four pairs at or below 0).
2. **Models versus ceiling** (26 pairs with ceiling above 0.05). Leave-publication-out R² is a median 55% of the ceiling and never exceeds it; scaffold-split R² exceeds the ceiling in 6 of 26 pairs, which a real model cannot do, consistent with leakage of publication offsets. Across pairs the ceiling predicts leave-publication-out R² (Spearman 0.73, p = 5e-6).
3. **Variance split.** About 57% of between-publication variance is protocol offset or noise; the chemical-series component is about 27% of total variance.
4. **Chronological split** (`05`). Train on publications up to the 75th-percentile year (median cut 2017), test on later ones, excluding test compounds from training: median R² 0.07 (45 pairs), within-series R² 0.02, lower than leave-publication-out.

5. **Empirical-Bayes anchor calibration** (`09`; the naive version in `07` failed, see below). Shrinking the few-anchor offset towards the training-set mean offset, with the between-publication variance estimated by *nested publication-level cross-validation*, beats the plain offset at every anchor count (pooled R², median over 47 pairs, plain to EB): 1 anchor 0.26 to 0.41 (better in 45/47 pairs, Wilcoxon p = 6e-13), 2 anchors 0.37 to 0.45 (42/47), 3 anchors 0.47 to 0.50 (36/47), 5 anchors 0.51 to 0.52 (33/47). One EB anchor therefore recovers more than a plain 2-anchor offset. Three splits, training capped at 2,500 compounds.
6. **Stronger models do not help** (`09`, 47 pairs, leave-publication-out, same capped training): pooled R² median RF 0.27, gradient boosting 0.17, MLP 0.17, 3-model mean 0.27; within-series R² 0.03, -0.13, -0.11 and 0.00; series-centred RF -0.03 pooled / 0.02 within. This reproduces the manuscript's conclusion (no learner improves within-series ranking).
7. **Equilibrium versus functional data** (`09`): median between-publication share 0.56 for Ki/Kd (37 pairs) and 0.66 for IC50 (45 pairs); both far above zero, so the publication effect is not an artefact of mixing assay types.
8. Independent validation of item 5 (`10`; four fresh splits, training capped at 6,000): plain vs EB pooled R² 0.18 vs 0.34 (1 anchor), 0.38 vs 0.46 (2), 0.44 vs 0.47 (3), 0.49 vs 0.50 (5); EB better in 46, 41, 41, 36 of 47 pairs. Figure: `results/figure_ceiling_and_calibration.png` (`11`).
9. Negative result: the naive out-of-bag version of shrinkage (`07`) underperforms the plain offset (0.28 vs 0.43 at 2 anchors) because out-of-bag residuals still carry same-publication information.

Per-pair tables: `results/*.csv`; headline numbers: `results/summary.json`, `results/model_summary.json`.

### Caveats
1. Repeat compounds are mostly reference or tool compounds and may not represent typical series compounds.
2. Vx lumps protocol offset with random error; they are not separated here.
3. Pairs with few repeats have wide CIs; summaries use pairs with at least 20 repeats.
4. The unfiltered pull gives 642k raw records before curation and 375k curated compound-document-target values; the manuscript reports 435,165 measurements for 237,568 compounds, and its exact filters are not recoverable, so counts are not expected to match.
5. The ceiling comparison uses repeat compounds, mostly reference compounds, and the scaffold analysis uses per-compound median selectivity, so scaffold R² and the ceiling are not strictly like-for-like.

## Reproduce
```
pip install -r requirements.txt
python code/01_fetch_chembl.py          # about 30 min, ~300 MB into raw/
python code/02_variance_components.py   # about 5 min
python code/03_summarise.py
python code/04_models_scaffold_vs_lpo.py raw results   # about 15 min on 2 CPUs
python code/05_time_split.py raw results
python code/06_ceiling_vs_models.py
python code/07_anchor_shrinkage.py raw results
python code/09_extra_models_eb_strat.py raw results
python code/10_eb_validation.py raw results
python code/11_figure.py
```
`data/curated_pX_by_compound_document_target.csv.gz` lets you skip the download: it holds the curated pX per compound, document and target (no SMILES; ChEMBL IDs only).

## Licence
Code: MIT. Derived data come from ChEMBL (CC BY-SA 3.0, EMBL-EBI); cite ChEMBL if you reuse them.
