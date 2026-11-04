# PARALOG: variance components of paralog-selectivity data

Code and derived data for the variance-component analysis that accompanies the manuscript on
leave-publication-out (LPO) evaluation of paralog/isoform selectivity models (ChEMBL 37, 51 target pairs).

## What this repository reproduces

The original 01–16 script pipeline lived in a sandbox that no longer exists, so the model analyses were **re-implemented** from the manuscript's Methods and re-run on a fresh ChEMBL 37 pull. The re-implementation matches the manuscript closely (below). Not re-implemented: XGBoost/MLP/GIN comparison, series-centred and pairwise training, Ki/Kd vs IC50 stratification, and the PDB pocket-divergence analysis.

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

5. **Naive shrinkage of the anchor offset does not help** (`07`). Shrinking the few-anchor offset by a variance ratio estimated from out-of-bag residuals gave lower pooled R² than the plain offset at 2, 3 and 5 anchors (median 0.28, 0.29, 0.27 vs 0.43, 0.48, 0.52; the plain offset was better in 39 to 47 of 47 pairs). Out-of-bag residuals still carry same-publication information, so the between-publication variance (and hence the shrinkage weight) is underestimated; a proper estimate needs nested publication-level cross-validation. Reported as a negative result.

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
```
`data/curated_pX_by_compound_document_target.csv.gz` lets you skip the download: it holds the curated pX per compound, document and target (no SMILES; ChEMBL IDs only).

## Licence
Code: MIT. Derived data come from ChEMBL (CC BY-SA 3.0, EMBL-EBI); cite ChEMBL if you reuse them.
