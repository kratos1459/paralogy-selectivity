"""Select an unbiased-by-hand set of target pairs: human single-protein ChEMBL targets that share the deepest ChEMBL protein-family class,
have >=MINREC exact IC50/Ki/Kd (nM) records, excluding the 51 original pairs; at most 3 random pairs per family class; 60 random pairs (seed 1).
Writes data_random/pairs.txt (keys are ChEMBL target IDs) and data_random/targets_resolved.json. Run from repo root."""
import requests,json,random,os,collections
from concurrent.futures import ThreadPoolExecutor
API='https://www.ebi.ac.uk/chembl/api/data/'; MINREC=500
def get(u,p=None):
    for i in range(5):
        try:
            r=requests.get(u,params=p,timeout=90)
            if r.status_code==200: return r.json()
        except Exception: pass
    raise RuntimeError(u)
pc={}; d=get(API+'protein_classification.json',{'limit':1000})
for x in d['protein_classifications']: pc[x['protein_class_id']]=x
print(len(pc),'classes',flush=True)
T=[]; d=get(API+'target.json',{'organism':'Homo sapiens','target_type':'SINGLE PROTEIN','limit':1000,'only':'target_chembl_id,pref_name,target_components'})
while True:
    T+=d['targets']
    if not d['page_meta']['next']: break
    d=get('https://www.ebi.ac.uk'+d['page_meta']['next'])
print(len(T),'targets',flush=True)
comp={}
for t in T:
    for c in t['target_components']: comp[c['component_id']]=None
cids=sorted(comp)
def cls(chunk):
    d=get(API+'target_component.json',{'component_id__in':','.join(map(str,chunk)),'limit':100,'only':'component_id,protein_classifications'})
    return {c['component_id']:[x['protein_classification_id'] for x in c['protein_classifications']] for c in d['target_components']}
with ThreadPoolExecutor(8) as ex:
    for r in ex.map(cls,[cids[i:i+100] for i in range(0,len(cids),100)]): comp.update(r)
groups=collections.defaultdict(list)
for t in T:
    cs=[x for c in t['target_components'] for x in (comp.get(c['component_id']) or [])]
    if cs: groups[max(cs,key=lambda k:pc[k]['class_level'])].append(t)
groups={k:v for k,v in groups.items() if len(v)>=2 and pc[k]['class_level']>=3}
print(len(groups),'candidate classes',flush=True)
def count(t):
    d=get(API+'activity.json',{'target_chembl_id':t['target_chembl_id'],'standard_type__in':'IC50,Ki,Kd','standard_units':'nM','standard_relation':'=','limit':1,'only':'activity_id'})
    return t['target_chembl_id'],d['page_meta']['total_count']
need=[t for v in groups.values() for t in v]
with ThreadPoolExecutor(8) as ex: n=dict(ex.map(count,need))
orig=set(frozenset(json.load(open('data/targets_resolved.json'))[g][0][0] for g in l.strip().split('/')) for l in open('data/pairs.txt'))
T0=json.load(open('data/targets_resolved.json'))
random.seed(1); cand=[]
for k,v in groups.items():
    v=[t for t in v if n.get(t['target_chembl_id'],0)>=MINREC]
    ps=[(a,b) for i,a in enumerate(v) for b in v[i+1:] if frozenset([a['target_chembl_id'],b['target_chembl_id']]) not in {frozenset([T0[g][0][0] for g in o]) if False else o for o in orig}]
    ps=[p for p in ps if frozenset([p[0]['target_chembl_id'],p[1]['target_chembl_id']]) not in orig]
    random.shuffle(ps); cand+=[(k,)+p for p in ps[:3]]
random.shuffle(cand); sel=cand[:60]
os.makedirs('data_random',exist_ok=True); tj={}
with open('data_random/pairs.txt','w') as f:
    for k,a,b in sel:
        ia,ib=a['target_chembl_id'],b['target_chembl_id']; tj[ia]=[[ia,a['pref_name']]]; tj[ib]=[[ib,b['pref_name']]]; f.write(f'{ia}/{ib}\n')
json.dump(tj,open('data_random/targets_resolved.json','w'),indent=1)
print(len(sel),'pairs;',len(tj),'targets; classes used',len({k for k,_,_ in sel}))
