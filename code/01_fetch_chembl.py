"""Resolve gene symbols (data/pairs.txt) to human ChEMBL SINGLE PROTEIN targets and download
IC50/Ki/Kd (nM) activities for biochemical (B) and functional (F) assays via the public ChEMBL web API.
Writes raw/<target>.jsonl (not committed: ~300 MB). Tested against ChEMBL_37 (2026-05-01).
Note: ADRA1A resolves to two targets; CHEMBL229 (alpha-1A) is selected explicitly."""
import requests,json,os,time
from concurrent.futures import ThreadPoolExecutor
API='https://www.ebi.ac.uk/chembl/api/data/'
genes=sorted({g for l in open('data/pairs.txt') for g in l.strip().split('/')})
T={}
for g in genes:
    r=requests.get(API+'target.json',params={'target_components__target_component_synonyms__component_synonym__iexact':g,
        'target_type':'SINGLE PROTEIN','organism':'Homo sapiens','limit':20},timeout=60).json()
    T[g]=[(t['target_chembl_id'],t['pref_name']) for t in r['targets']]
T['ADRA1A']=[('CHEMBL229','Alpha-1A adrenergic receptor')]
json.dump(T,open('data/targets_resolved.json','w'),indent=1)
tids=sorted({v[0][0] for v in T.values()})
os.makedirs('raw',exist_ok=True)
F=('molecule_chembl_id,canonical_smiles,standard_type,standard_value,standard_units,standard_relation,document_chembl_id,'
   'assay_chembl_id,pchembl_value,data_validity_comment,potential_duplicate,target_chembl_id,assay_type')
def get(u,params=None):
    for i in range(6):
        try:
            r=requests.get(u,params=params,timeout=120)
            if r.status_code==200: return r.json()
        except Exception: pass
        time.sleep(3*(i+1))
    raise RuntimeError(u)
def run(t):
    fn=f'raw/{t}.jsonl'
    if os.path.exists(fn+'.done'): return
    n=0
    with open(fn,'w') as f:
        d=get(API+'activity.json',{'target_chembl_id':t,'standard_type__in':'IC50,Ki,Kd','standard_units':'nM',
              'assay_type__in':'B,F','only':F,'limit':1000})
        while True:
            for a in d['activities']: f.write(json.dumps(a)+'\n'); n+=1
            if not d['page_meta']['next']: break
            d=get('https://www.ebi.ac.uk'+d['page_meta']['next'])
    open(fn+'.done','w').write(str(n)); print(t,n,flush=True)
with ThreadPoolExecutor(8) as ex: list(ex.map(run,tids))
