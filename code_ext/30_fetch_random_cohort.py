#!/usr/bin/env python3
"""Stage 30: fetch the 60 randomly sampled ChEMBL target pairs and build the same derived tables as
the curated cohort, so the extension analyses can be replicated on them.

Pair list and target identifiers are those of the released random cohort (data_random/pairs.txt of
the repository). Activities are retrieved with the same query as code/01_fetch_chembl.py:
IC50/Ki/Kd in nM, all assay types. Writes data/raw_random/<target>.jsonl, then
data/doc_level_random.csv.gz, data/smiles_random.csv.gz and data/pairs_random.csv.

usage: python 30_fetch_random_cohort.py
"""
import json, os, time, requests
import numpy as np, pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw_random"; RAW.mkdir(parents=True, exist_ok=True)
API = "https://www.ebi.ac.uk/chembl/api/data/"
PAIRS = """CHEMBL3892/CHEMBL2955 CHEMBL211/CHEMBL1821 CHEMBL4527/CHEMBL6166 CHEMBL4895/CHEMBL5314 CHEMBL249/CHEMBL4429
CHEMBL4478/CHEMBL4641 CHEMBL2028/CHEMBL1917 CHEMBL2148/CHEMBL2835 CHEMBL4142/CHEMBL2742 CHEMBL2069/CHEMBL1836
CHEMBL2250/CHEMBL1841 CHEMBL209/CHEMBL244 CHEMBL1741207/CHEMBL1741215 CHEMBL3706/CHEMBL2285 CHEMBL6032/CHEMBL2189110
CHEMBL4618/CHEMBL2693 CHEMBL4462/CHEMBL4506 CHEMBL2815/CHEMBL4898 CHEMBL226/CHEMBL256 CHEMBL2327/CHEMBL4429
CHEMBL1937/CHEMBL1829 CHEMBL233/CHEMBL2014 CHEMBL224/CHEMBL3371 CHEMBL4662/CHEMBL5620 CHEMBL1785/CHEMBL252
CHEMBL1829/CHEMBL325 CHEMBL2147/CHEMBL5407 CHEMBL1941/CHEMBL231 CHEMBL4852/CHEMBL5330 CHEMBL1841/CHEMBL3905
CHEMBL4179/CHEMBL2637 CHEMBL219/CHEMBL234 CHEMBL268/CHEMBL3837 CHEMBL3807/CHEMBL2889 CHEMBL4441/CHEMBL2010631
CHEMBL210/CHEMBL223 CHEMBL1898/CHEMBL225 CHEMBL227/CHEMBL4607 CHEMBL1913/CHEMBL1936 CHEMBL4608/CHEMBL3795
CHEMBL2996/CHEMBL3920 CHEMBL4625/CHEMBL4860 CHEMBL1741195/CHEMBL2176771 CHEMBL4393/CHEMBL333 CHEMBL1941/CHEMBL264
CHEMBL3227/CHEMBL5137 CHEMBL2335/CHEMBL3976 CHEMBL2434/CHEMBL2010631 CHEMBL1868/CHEMBL1955 CHEMBL2107/CHEMBL2010631
CHEMBL1881/CHEMBL3710 CHEMBL2889/CHEMBL4903 CHEMBL1795087/CHEMBL2157850 CHEMBL4895/CHEMBL5331 CHEMBL4015/CHEMBL4423
CHEMBL4393/CHEMBL4073 CHEMBL217/CHEMBL219 CHEMBL4801/CHEMBL3468 CHEMBL4816/CHEMBL4282 CHEMBL4439/CHEMBL5311""".split()
F = ("molecule_chembl_id,canonical_smiles,standard_type,standard_value,standard_units,standard_relation,"
     "document_chembl_id,assay_chembl_id,pchembl_value,data_validity_comment,potential_duplicate,"
     "target_chembl_id,assay_type")


def get(u, params=None):
    for i in range(6):
        try:
            r = requests.get(u, params=params, timeout=120)
            if r.status_code == 200: return r.json()
        except Exception:
            pass
        time.sleep(3 * (i + 1))
    raise RuntimeError(u)


def run(t):
    fn = RAW / f"{t}.jsonl"
    if (RAW / f"{t}.jsonl.done").exists(): return
    n = 0
    with open(fn, "w") as f:
        d = get(API + "activity.json", {"target_chembl_id": t, "standard_type__in": "IC50,Ki,Kd",
                                        "standard_units": "nM", "only": F, "limit": 1000})
        while True:
            for a in d["activities"]:
                f.write(json.dumps(a) + "\n"); n += 1
            if not d["page_meta"]["next"]: break
            d = get("https://www.ebi.ac.uk" + d["page_meta"]["next"])
    open(str(fn) + ".done", "w").write(str(n)); print(t, n, flush=True)


tids = sorted({x for p in PAIRS for x in p.split("/")})
print(len(tids), "random-cohort targets")
with ThreadPoolExecutor(8) as ex: list(ex.map(run, tids))

rows, smi = [], {}
for t in tids:
    fn = RAW / f"{t}.jsonl"
    if not fn.exists(): continue
    df = pd.read_json(fn, lines=True)
    if df.empty: continue
    df = df[(df.standard_relation == "=") & df.standard_value.notna() & df.canonical_smiles.notna()]
    df = df[df.data_validity_comment.isna() & (df.potential_duplicate.fillna(0).astype(int) == 0)]
    if df.empty: continue
    df["pX"] = 9 - np.log10(df.standard_value.astype(float).clip(lower=1e-6))
    df = df[(df.pX >= 2) & (df.pX <= 12.5)]
    smi.update(df.drop_duplicates("molecule_chembl_id").set_index("molecule_chembl_id").canonical_smiles.to_dict())
    r = df.groupby(["molecule_chembl_id", "document_chembl_id"]).pX.median().reset_index()
    r["target_chembl_id"] = t
    rows.append(r[["target_chembl_id", "molecule_chembl_id", "document_chembl_id", "pX"]])

D = pd.concat(rows, ignore_index=True)
D.to_csv(ROOT / "data" / "doc_level_random.csv.gz", index=False)
pd.DataFrame({"molecule_chembl_id": list(smi), "canonical_smiles": list(smi.values())}).to_csv(
    ROOT / "data" / "smiles_random.csv.gz", index=False)
P = pd.DataFrame([dict(family="random", A=a, tidA=a, B=b, tidB=b) for a, b in (p.split("/") for p in PAIRS)])
P.to_csv(ROOT / "data" / "pairs_random.csv", index=False)
print("doc_level_random", D.shape, "smiles", len(smi), "pairs", len(P))
