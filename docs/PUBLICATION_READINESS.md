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
* Empirical-Bayes anchor calibration with nested publication-level CV: significantly better than the plain offset at every anchor count (1 anchor: R² 0.26 to 0.41). This is the one method-level result and the strongest candidate for a higher tier; it should be validated on held-out pairs/families and with a fully uncapped training set.
* Stronger learners (GBM, MLP, ensemble, series-centred) do not beat RF; Ki/Kd vs IC50 stratification keeps the publication effect (0.56 / 0.66).
* Tried to separate random error from protocol offset using within-publication record variance: inconclusive (repeated records within a publication mix different assay variants), so not shipped. Needs assay descriptors.

## Still requires something other than conversation-level compute
* Prospective or single-laboratory data (wet lab or an external panel), pocket-divergence and stronger-model analyses re-run in the repo (possible with more compute time).

## Known weaknesses a reviewer will raise
* Documents are a proxy for series and laboratories.
* Repeat compounds are not representative of series compounds.
* Pairs are curated, not random, and share targets.
* Re-implemented model results differ from the original by up to 0.05 in R².

## What is still needed for a high-impact journal (checklist)
1. **External validity (the main gap).** One of: (a) single-laboratory or standardised-panel data (kinome or safety-panel profiling) showing the ceiling and calibration gains; (b) prospective measurement of a few anchors in a collaborator's lab; (c) a curated external set from after the ChEMBL 37 cutoff.
2. **Method claim with theory.** Formalise the empirical-Bayes anchor offset (hierarchical model, shrinkage weight, expected R2 gain) and test on pairs and families not used to tune anything.
3. **Complete, reproducible pipeline.** Re-run XGBoost, a graph network and a pretrained model; pocket-divergence analysis; use the same pipeline for every number in the paper.
4. **Benchmark release.** Splits, curation code, leaderboard baselines, a Zenodo DOI and a data paper style description.
5. **Separate protocol offset from error.** Needs assay descriptors (cell line, enzyme construct, substrate, assay format) from ChEMBL assay records.
6. **Positioning and novelty check.** Re-run the literature search before submission; position against matched-pair, assay-noise and temporal-split work.
7. **Reporting.** Drop the "chosen pairs" limitation by adding a random sample of pairs, add per-family results and confidence intervals everywhere.
Realistic targets: J Cheminform or JCIM now; Briefings in Bioinformatics or Nature Communications only after items 1 to 4.

## Update after the three extra analyses
* Theory: closed-form shrinkage predictions match observed error (r 0.95 to 0.99), which turns the calibration from a heuristic into a model-based result.
* Assay descriptors: recorded descriptors explain only about 20% of cross-document disagreement, so most of the noise ceiling is unrecorded laboratory effects (supports the need for true single-lab or prospective data).
* Single-document hold-out: base models fail (median R2 -0.11); EB calibration still helps in 42 to 44 of 47 tests. This is the honest weak spot; a reviewer will read it as the external-validity gap. The paper should state it openly and frame calibration as the practical remedy.
* Still not possible here: genuinely single-laboratory or prospective data, and a Zenodo DOI (needs your account).

## Update (7 Nov)
* XGBoost re-run (real xgboost): worse than RF (0.14 vs 0.20 pooled; within-series -0.14), matching the manuscript.
* Benchmark splits released in `benchmark/`; theory derivation in `docs/THEORY.md`.
* Open items that need you: Zenodo DOI (link the GitHub repo in Zenodo), true single-lab or prospective data, graph-network and pocket-divergence re-runs.
