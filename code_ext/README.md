# Extension analyses (stages 25-29)

Added after the ACS Omega desk rejection, which asked for the chemistry behind the evaluation
result. Each script reads the released derived data (`data/`) and writes to `results_ext/` and
`figures_ext/`; `common_ext.py` rebuilds the matched same-publication endpoint and reproduces the
leave-publication-out protocol of `code/04_models_scaffold_vs_lpo.py`.

| script | question | key output |
|---|---|---|
| `25_representations.py` | Does the molecular representation change cross-publication generalisation? Morgan r2/r3, atom pairs, torsions, MACCS, 208 RDKit descriptors, Morgan+descriptors. | `representation_benchmark_*.csv` |
| `26_matched_molecular_pairs.py` | Do single substituent changes carry reproducible selectivity effects between laboratories? Within-publication MMPs, cross-publication replication stratified by effect size, transfer model. | `mmp_*.csv`, `mmp_summary.json` |
| `27_physchem_drivers.py` | Which physicochemical properties drive selectivity once the publication offset is removed? Twelve properties, publication-centred, permutation p-values. | `physchem_*.csv` |
| `28_pocket_substitution_chemistry.py` | Which pocket substitutions open a selectivity window? Volume, hydropathy, charge, aromaticity and H-bonding changes at differing pocket positions. | `pocket_chemistry_*.csv`, `pocket_substitutions.csv` |
| `29_figures_ext.py` | Figures E1-E3 at 600 dpi (PNG and LZW TIFF). | `figures_ext/` |

Headline results and ready-to-paste manuscript text are in `EXTENSION_TEXT.md`.

```bash
cd code_ext
python 26_matched_molecular_pairs.py
python 27_physchem_drivers.py
python 28_pocket_substitution_chemistry.py
python 25_representations.py            # or shard: python 25_representations.py 0 17 _s0
python 29_figures_ext.py
```
