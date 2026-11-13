# Extension analyses (stages 25-29): the chemistry behind the evaluation result

Written after the ACS Omega desk rejection, whose stated reason was that the work "fails to relate
its results with the chemistry behind them". These four analyses answer that objection with
chemistry rather than statistics: how large a substituent effect has to be before it survives
between laboratories, what physicochemistry drives paralog selectivity once the publication offset
is removed, and which pocket substitutions open a selectivity window.

## Methods — new paragraph

Matched molecular pairs, physicochemical drivers, pocket substitution chemistry, and molecular
representations. Each molecule was cut at every acyclic single bond between heavy atoms, giving
(core, R-group) decompositions; two molecules measured in the same publication that share a core and
differ in one R-group form a matched molecular pair (MMP) whose effect is the difference of their
same-publication selectivities, so the publication offset cancels by construction. Effects were
aggregated per (transformation, publication), and transformations observed in at least two
publications were used to measure cross-publication reproducibility (Spearman correlation, sign
concordance and the ratio of replicate effect magnitudes), stratified by the magnitude of the effect
in the first publication. A transfer model predicted each MMP effect in a held-out publication from
the median effect of the same transformation in the remaining publications and was scored against
predicting no effect. Physicochemical drivers were assessed for twelve RDKit properties by centring
both selectivity and property within each publication and computing Spearman correlations with
publication-wise permutation p-values (500 permutations), together with the within-series variance
explained by the best single property and by all properties jointly. Binding-pocket substitutions
were derived from the cached pocket residue sets and UniProt sequences of the main analysis: every
pocket position differing between the two paralogs was classified by its change in side-chain volume,
Kyte-Doolittle hydropathy, formal charge, aromaticity and hydrogen-bonding capacity, and the
resulting descriptors were related to the selectivity window by Spearman correlation with 2,000-fold
bootstrap confidence intervals and 5,000 within-family permutations. Finally, seven molecular
representations (Morgan counts of radius 2 and 3, atom pairs, topological torsions, MACCS keys, 208
RDKit physicochemical descriptors, and Morgan concatenated with those descriptors) were run through
the identical leave-publication-out protocol with the same random forest.

## Results — new subsection: "Substituent effects below ten-fold do not survive between laboratories"

Matched molecular pairs remove the publication offset almost completely: across 46 target pairs with
sufficient matched pairs, the between-publication share of variance fell from a median of 0.65 for
raw selectivity to 0.11 for transformation effects, lower in 46 of 46 pairs (paired Wilcoxon
p = 3 x 10^-14; Figure E1a). What remains, however, is barely reproducible. Across 77,417 instances
in which the same transformation was measured in two independent publications, effects correlated
only weakly (pooled Spearman 0.34; median per-pair 0.14) and agreed in sign for 62% of instances
(median per pair 54%). A transfer model that predicts a held-out publication's MMP effects from the
median effect of the same transformation elsewhere was worse than predicting no effect at all in 42
of 46 pairs (median R^2 -0.49).

Stratifying by effect size turns this into a usable rule. Transformations whose effect in the first
publication was below 0.3 log units reproduced in sign only 54% of the time - indistinguishable from
chance - and their apparent magnitude grew on replication (median |effect| 0.13 to 0.23, ratio 1.80),
the signature of regression to noise. Sign concordance rose to 65% at 0.3-0.5 log units, 69% at
0.5-1.0 and 71% at 1.0-2.0, while the magnitude shrank steadily (ratios 0.85, 0.66 and 0.52;
Figure E1c). A median transformation shifts selectivity by 0.60 log units, against a cross-publication
replicate noise of comparable size, so most single-point modifications are simply not resolvable
between laboratories. Only changes of roughly one log unit or more - about ten-fold - carry a
reproducible direction, and even those arrive at about half their originally reported magnitude.

## Results — new subsection: "Selectivity is pocket-local, not a bulk property"

Once selectivity and properties are centred within each publication, bulk physicochemistry explains
almost none of the within-series signal: the median absolute correlation was 0.054 or less for every
one of twelve properties (molecular weight, cLogP, TPSA, hydrogen-bond donors and acceptors,
rotatable bonds, aromatic and aliphatic rings, fraction sp3, heavy atoms, formal charge, molar
refractivity), against 0.13-0.17 for the same correlations computed without centring - so the
apparent property dependence of selectivity is largely publication composition. The best single
property explained a median of 1.7% of within-series selectivity variance, and all twelve jointly
4.8%; no property dominated across pairs (fraction sp3 was best in 9 of 47 pairs, TPSA in 6,
hydrogen-bond donors in 5, cLogP in 5; Figure E2a,b).

The selectivity window is instead set by how many pocket positions differ. Across 40 pairs with
pocket structures for both targets, the spread of same-publication selectivity correlated with the
number of differing pocket positions (Spearman 0.62, within-family permutation p < 0.001), with the
summed side-chain volume change (0.58), the summed hydropathy change (0.52), the number of
hydrogen-bonding changes (0.49) and the number of aromaticity changes (0.46), but not with the mean
severity of a substitution (mean |hydropathy change| -0.04, p = 0.83; Figure E2c). Selectivity
windows are therefore opened by the accumulation of many modest pocket differences rather than by
one drastic substitution, which is consistent with the weak, distributed structure-selectivity signal
seen in the matched-pair analysis.

## Discussion — paragraph to add

These analyses locate the problem in the chemistry rather than in the statistics. A median
single-substituent change moves paralog selectivity by about 0.6 log units, which is the same size as
the disagreement between two laboratories measuring the same compound pair; below roughly one log
unit, the direction of a reported substituent effect is no better than a coin toss on replication,
and reported magnitudes shrink by half when they are reproduced. That is why no amount of model
capacity or feature engineering rescues cross-publication prediction: the within-series signal that a
model would have to learn is itself mostly irreproducible. It also explains the pocket result -
selectivity windows widen with the number of differing pocket positions rather than with the severity
of any single substitution - and the failure of bulk physicochemical descriptors, since selectivity
arises from local complementarity rather than from global lipophilicity or size. For practice, the
implication is concrete: treat reported selectivity shifts under about ten-fold as unconfirmed until
they are reproduced in one laboratory, and design selectivity campaigns around pocket positions that
differ in number, not around a single hoped-for contact.

## Figures

Figure E1. Matched molecular pair analysis. (a) Between-publication share of variance for raw
selectivity against the same share for MMP transformation effects, one point per target pair, coloured
by family; the dashed line is equality. (b) Effect of the same transformation measured in two
independent publications (one point per transformation; axes clipped at +/-3 log units). (c) Sign
concordance between publications (bars, left axis) and the ratio of replicate effect magnitudes (line,
right axis) by the size of the effect in the first publication.

Figure E2. Chemistry of selectivity. (a) Median absolute Spearman correlation between selectivity and
each physicochemical property, computed raw and within publication. (b) Distribution across pairs of
the within-series variance explained by the best single property. (c) Spearman correlation of each
pocket substitution descriptor with the selectivity window (bars; error bars are 2,000-fold bootstrap
confidence intervals; red bars are not significant under within-family permutation).

Figure E3. Representation benchmark. Leave-publication-out and within-series R^2 for seven molecular
representations under the identical protocol.

## Suggested reframing for resubmission

The present title leads with the evaluation protocol, which is what both editors read as "not
chemistry". A chemistry-led framing of the same work, with the new analyses as the opening result:

  Title: "How Large Must a Selectivity Change Be to Reproduce? Matched-Pair and Binding-Pocket
  Analysis of Paralog Selectivity Across 51 Human Target Pairs"

  Abstract first sentence: "A median single-substituent change shifts selectivity between two human
  paralogs by 0.6 log units, the same size as the disagreement between two laboratories measuring the
  same compounds, so only changes above roughly ten-fold reproduce in direction and even those arrive
  at half their reported magnitude."

The evaluation-protocol results (scaffold versus leave-publication-out, the noise ceiling, the
empirical-Bayes anchor estimator) then follow as the explanation of why models trained on such data do
not transfer, instead of being the headline.

## Numbers to cite

| quantity | value |
|---|---|
| target pairs with matched molecular pairs | 46 |
| between-publication variance share: selectivity vs MMP effect | 0.65 vs 0.11 (46/46 pairs, p = 3e-14) |
| cross-publication replicate instances | 77,417 |
| pooled Spearman between replicate effects | 0.34 (median per pair 0.14) |
| pooled sign concordance | 62% (median per pair 54%) |
| MMP transfer model | median R^2 -0.49; better than predicting zero in 4/46 pairs |
| sign concordance by effect size | 54% (<0.3), 65% (0.3-0.5), 69% (0.5-1.0), 71% (1.0-2.0) |
| magnitude ratio on replication | 1.80, 0.85, 0.66, 0.52, 0.23 across the same strata |
| median SD of a transformation effect / of selectivity | 0.60 / 1.04 log units |
| within-publication property correlations | median |rho| <= 0.054 for all 12 properties (raw 0.13-0.17) |
| within-series variance explained by properties | 1.7% best single, 4.8% all twelve |
| pocket descriptors vs selectivity window | n differing positions 0.62, sum |dVolume| 0.58, sum |dHydropathy| 0.52, H-bond changes 0.49, aromatic changes 0.46, mean |dHydropathy| -0.04 (ns) |
