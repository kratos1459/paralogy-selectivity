# PARALOG: variance components of paralog-selectivity data

Code and derived data for the variance-component analysis that accompanies the manuscript on
leave-publication-out (LPO) evaluation of paralog/isoform selectivity models (ChEMBL 37, 51 target pairs).

## What this repository reproduces

| Item | Status |
|---|---|
| Pair list (51 pairs, from SI Table S1) | `data/pairs.txt`; gene to ChEMBL target mapping in `data/targets_resolved.json` |
| ChEMBL download | `code/01_fetch_chembl.py` (public API, ChEMBL_37, 2026-05-01) |
| Between-publication share of selectivity variance | reproduced (`code/02_variance_components.py`) |
| Protocol/replicate noise and noise ceiling | new analysis, same script |
| Model results (scaffold-split R2 0.56, LPO R2 0.21, within-series 0.05, anchor calibration) | **not included**: the original model-training pipeline is not part of this release |

The earlier files `selectivity.csv` and `features.npz` were not deposited; they contained non-structure tokens in the SMILES column.

## Results (50 of 51 pairs analysable; CACNA1C/CACNA1D has no same-publication overlap)

* Between-publication share of delta variance: median 0.69 (IQR 0.59 to 0.78); manuscript reports 0.65 (0.54 to 0.76).
* Cross-publication replicate variance of the same compound's delta (protocol plus measurement noise), 30 pairs with at least 20 repeat compounds: median 41% of total variance (IQR 25% to 64%).
* Noise ceiling for R2 on a held-out publication, 1 - Vx/Vtot: median 0.59 (IQR 0.36 to 0.75); four pairs have a ceiling at or below 0 (BRD4/BRD2, EGFR/ERBB2, CDK2/CDK1, CA2/CA1).
* Of the between-publication variance, a median 58% is explained by protocol offsets; the remainder (median 27% of total variance) is chemical-series effect.
* Per-pair values and bootstrap 95% CIs: `results/varcomp_by_pair.csv`; summary: `results/summary.json`.

### Caveats
1. Repeat compounds are mostly reference or tool compounds and may not represent typical series compounds.
2. Vx lumps protocol offset with random error; they are not separated here.
3. Pairs with few repeats have wide CIs; summaries use pairs with at least 20 repeats.
4. The raw pull (635,280 activity records before curation) gives 375,452 curated compound-document-target values; this differs from the 435,165 measurements reported in the manuscript, whose exact filters were not available, so counts are not expected to match.

## Reproduce
```
pip install -r requirements.txt
python code/01_fetch_chembl.py          # about 30 min, ~300 MB into raw/
python code/02_variance_components.py   # about 5 min
python code/03_summarise.py
```
`data/curated_pX_by_compound_document_target.csv.gz` lets you skip the download: it holds the curated pX per compound, document and target (no SMILES; ChEMBL IDs only).

## Licence
Code: MIT. Derived data come from ChEMBL (CC BY-SA 3.0, EMBL-EBI); cite ChEMBL if you reuse them.
