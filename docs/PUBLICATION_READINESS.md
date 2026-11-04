# Publication readiness (candid assessment)

## What the data now support
* Scaffold-split performance overstates generalisation; the effect reproduces independently (0.54 vs 0.26 median R²; 92% leakage; 65% between-publication variance).
* New: a measurable **noise ceiling** (median R² 0.63). Publication-held-out models reach about 55% of it and never exceed it; scaffold-split models exceed it in 6 of 26 pairs, which is evidence of leakage rather than skill.
* New: ceiling predicts learnability across pairs (Spearman 0.73) and a **chronological split** is harsher still (median R² 0.07).

## Realistic venue
* J Cheminform / JCIM / J Med Chem (benchmarking): plausible with the current evidence.
* Higher tier (e.g. Nature Communications, Briefings in Bioinformatics, Nature Machine Intelligence): would need all of the following.
  1. A single-laboratory or standardised-panel validation (e.g. kinome profiling sets) showing the ceiling is protocol-driven rather than a ChEMBL curation artefact.
  2. Separate protocol offset from random error (needs records with assay descriptors, or replicate designs).
  3. A released benchmark with the publication-held-out and chronological splits, leaderboard-style baselines (tuned GBM, pretrained or pairwise models), and a citable DOI.
  4. A method result: for example empirical-Bayes anchor shrinkage that beats offset calibration under the same hold-out.
  5. Re-run of the remaining manuscript analyses (stronger models, stratification, pocket divergence) in this repository so every number is reproducible.

## Done computationally since the first draft
* Re-implemented scaffold/LPO/anchor pipeline (matches manuscript), noise ceiling, chronological split.
* Tried empirical-Bayes anchor shrinkage: negative (naive version); nested-CV version not yet done.
* Tried to separate random error from protocol offset using within-publication record variance: inconclusive (repeated records within a publication mix different assay variants), so not shipped. Needs assay descriptors.

## Still requires something other than conversation-level compute
* Prospective or single-laboratory data (wet lab or an external panel), pocket-divergence and stronger-model analyses re-run in the repo (possible with more compute time).

## Known weaknesses a reviewer will raise
* Documents are a proxy for series and laboratories.
* Repeat compounds are not representative of series compounds.
* Pairs are curated, not random, and share targets.
* Re-implemented model results differ from the original by up to 0.05 in R².
