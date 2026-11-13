"""Shared helpers for the extension analyses (stages 25+).

Rebuilds the matched same-publication selectivity endpoint from the released derived data
(data/doc_level.csv.gz, data/smiles.csv.gz, data/pairs.csv) and reproduces the leave-publication-out
protocol of code/04_models_scaffold_vs_lpo.py: documents with at least 5 matched compounds, pairs
with at least 8 such documents and at least 60 matched compounds, 25% of documents held out, and
every molecule appearing in a held-out document purged from training.
"""
import numpy as np, pandas as pd
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator as G
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results_ext"
FIG = ROOT / "figures_ext"
RES.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
MIN_DOC_COMPOUNDS, MIN_DOCS, MIN_MOLS = 5, 8, 60


def load_all():
    doc = pd.read_csv(DATA / "doc_level.csv.gz")
    smi = pd.read_csv(DATA / "smiles.csv.gz").set_index("molecule_chembl_id").canonical_smiles.to_dict()
    pairs = pd.read_csv(DATA / "pairs.csv")
    return doc, smi, pairs


def matched(doc, tA, tB):
    a = doc[doc.target_chembl_id == tA][["molecule_chembl_id", "document_chembl_id", "pX"]]
    b = doc[doc.target_chembl_id == tB][["molecule_chembl_id", "document_chembl_id", "pX"]]
    m = a.merge(b, on=["molecule_chembl_id", "document_chembl_id"], suffixes=("A", "B"))
    if m.empty: return m
    return m.assign(d=m.pXA - m.pXB)


def lpo_frame(m):
    if m is None or len(m) == 0: return None
    dc = m.groupby("document_chembl_id").size()
    good = dc[dc >= MIN_DOC_COMPOUNDS].index
    mg = m[m.document_chembl_id.isin(good)]
    if len(good) < MIN_DOCS or mg.molecule_chembl_id.nunique() < MIN_MOLS: return None
    return mg.reset_index(drop=True)


def lpo_splits(mg, n_splits=3, seed0=100, frac=0.25):
    dd = mg.document_chembl_id.values; mol = mg.molecule_chembl_id.values
    udocs = np.array(sorted(set(dd)))
    for s in range(n_splits):
        rng = np.random.default_rng(seed0 + s)
        te_d = set(rng.choice(udocs, max(1, int(frac * len(udocs))), replace=False))
        te = np.isin(dd, list(te_d))
        tr = (~te) & (~np.isin(mol, list(set(mol[te]))))
        if te.sum() < 20 or tr.sum() < 30: continue
        yield s, tr, te


def centered(v, groups):
    return pd.Series(np.asarray(v)).groupby(np.asarray(groups)).transform(lambda s: s - s.mean()).values


def metrics(y, p, docs):
    y = np.asarray(y, float); p = np.asarray(p, float)
    yc, pc = centered(y, docs), centered(p, docs)
    sy, sp = yc.std(), pc.std()
    k = sp / sy if sy > 0 else np.nan
    r = float(np.corrcoef(yc, pc)[0, 1]) if sy > 0 and sp > 0 else np.nan
    return dict(r2=float(r2_score(y, p)), r2_within=float(r2_score(yc, pc)) if sy > 0 else np.nan,
                rho_within=float(spearmanr(yc, pc)[0]) if sy > 0 else np.nan, k=float(k), r_centred=r)


_mg2 = G.GetMorganGenerator(radius=2, fpSize=2048)
_mg3 = G.GetMorganGenerator(radius=3, fpSize=4096)
_ap = G.GetAtomPairGenerator(fpSize=2048)
_tt = G.GetTopologicalTorsionGenerator(fpSize=2048)
_molcache = {}
_featcache = {}


def _mol(mid, smi):
    if mid not in _molcache:
        s = smi.get(mid)
        _molcache[mid] = Chem.MolFromSmiles(s) if isinstance(s, str) and s else None
    return _molcache[mid]


def _desc(mol):
    from rdkit.Chem import Descriptors
    global _DFNS
    try: fns = _DFNS
    except NameError:
        fns = _DFNS = [f for _, f in Descriptors._descList]
    v = np.array([f(mol) for f in fns], dtype=np.float64)
    v = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)
    return np.clip(v, -1e6, 1e6).astype(np.float32)


def _one(mid, smi, rep):
    key = (mid, rep)
    if key in _featcache: return _featcache[key]
    from rdkit.Chem import MACCSkeys
    mol = _mol(mid, smi); v = None
    if mol is not None:
        try:
            if rep == "morgan2_2048": v = _mg2.GetCountFingerprintAsNumPy(mol).astype(np.float32)
            elif rep == "morgan3_4096": v = _mg3.GetCountFingerprintAsNumPy(mol).astype(np.float32)
            elif rep == "atompair_2048": v = _ap.GetCountFingerprintAsNumPy(mol).astype(np.float32)
            elif rep == "torsion_2048": v = _tt.GetCountFingerprintAsNumPy(mol).astype(np.float32)
            elif rep == "maccs_167": v = np.array(MACCSkeys.GenMACCSKeys(mol), dtype=np.float32)
            elif rep == "rdkit_descriptors": v = _desc(mol)
            elif rep == "morgan2_plus_descriptors":
                v = np.concatenate([_mg2.GetCountFingerprintAsNumPy(mol).astype(np.float32), _desc(mol)])
            else: raise ValueError(rep)
        except Exception:
            v = None
    _featcache[key] = v
    return v


def featurise(mids, smi, rep):
    """Return (matrix, keep_mask) aligned with mids."""
    keep = np.array([_one(m, smi, rep) is not None for m in mids])
    X = np.stack([_one(m, smi, rep) for m, k in zip(mids, keep) if k]) if keep.any() else np.zeros((0, 1), np.float32)
    return X, keep


REPS = ["morgan2_2048", "morgan3_4096", "atompair_2048", "torsion_2048", "maccs_167",
        "rdkit_descriptors", "morgan2_plus_descriptors"]
