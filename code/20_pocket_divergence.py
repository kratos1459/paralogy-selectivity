"""Binding-pocket divergence between paralog pairs (PDB ligand-bound structures) vs selectivity spread and LPO R2.
Usage: python 20_pocket_divergence.py OUTDIR  (needs network: ChEMBL, RCSB, PDBe/UniProt). Caches in OUTDIR/cache."""
import sys,os,json,time,re,numpy as np,pandas as pd,requests
from concurrent.futures import ThreadPoolExecutor
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Align import PairwiseAligner,substitution_matrices
from scipy.stats import spearmanr,rankdata
OUT=sys.argv[1]; os.makedirs(OUT+'/cache',exist_ok=True)
T=json.load(open('data/targets_resolved.json')); T['ADRA1A']=[['CHEMBL229','a']]
def get(u,**kw):
    for i in range(6):
        try:
            r=requests.get(u,timeout=60,**kw)
            if r.status_code==200: return r
        except Exception: pass
        time.sleep(2*(i+1))
    return None
def cached(fn,key,f):
    p=f'{OUT}/cache/{key}'
    if os.path.exists(p): return json.load(open(p))
    v=f(); 
    if v is not None: json.dump(v,open(p,'w'))
    return v
def accession(g):
    t=T[g][0][0]
    def f():
        r=get(f'https://www.ebi.ac.uk/chembl/api/data/target/{t}.json')
        return r.json()['target_components'][0]['accession'] if r else None
    return cached(f,f'acc_{g}.json',f)
def seq(acc):
    def f():
        r=get(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta') or get(f'https://www.ebi.ac.uk/proteins/api/proteins/{acc}',headers={'Accept':'text/x-fasta'})
        return ''.join(r.text.split('\n')[1:]) if r else None
    return cached(f,f'seq_{acc}.json',f)
def search(acc,res,methods):
    q={"query":{"type":"group","logical_operator":"and","nodes":[
      {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession","operator":"exact_match","value":acc}},
      {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.nonpolymer_entity_count","operator":"greater","value":0}},
      {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","operator":"less_or_equal","value":res}},
      {"type":"terminal","service":"text","parameters":{"attribute":"exptl.method","operator":"in","value":methods}}]},
      "return_type":"entry","request_options":{"paginate":{"start":0,"rows":60},"sort":[{"sort_by":"rcsb_entry_info.resolution_combined","direction":"asc"}]}}
    for i in range(5):
        try:
            r=requests.post('https://search.rcsb.org/rcsbsearch/v2/query',json=q,timeout=60)
            if r.status_code==200: return [x['identifier'] for x in r.json().get('result_set',[])]
            if r.status_code==204: return []
        except Exception: pass
        time.sleep(2*(i+1))
    return []
EXCL=set('GOL EDO PEG PGE PG4 1PE P6G SO4 PO4 ACT ACY DMS FMT TRS MPD BME CL NA K MG ZN CA MN FE CU CD NI CO IOD BR EPE MES HEZ CIT TAR NH4 NO3 SCN AZI IMD ATP ADP AMP ANP GTP GDP GNP NAD NAP NDP FAD FMN HEM HEC SAM SAH COA PLP NAG BMA MAN FUC GAL GLC NDG BGC SIA UDP CLR OLC LMT BOG DAO HTO PLM OLA STE MYR'.split())
AA=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
def pocket(pdb,acc):
    def f():
        r=get(f'https://files.rcsb.org/download/{pdb}.cif')
        if not r: return None
        open(f'/tmp/{pdb}.cif','w').write(r.text); d=MMCIF2Dict(f'/tmp/{pdb}.cif'); os.remove(f'/tmp/{pdb}.cif')
        def L(k): v=d.get(k,[]); return v if isinstance(v,list) else [v]
        refs={}
        for a,s,b,ab in zip(L('_struct_ref_seq.pdbx_db_accession'),L('_struct_ref_seq.pdbx_strand_id'),L('_struct_ref_seq.db_align_beg'),L('_struct_ref_seq.pdbx_auth_seq_align_beg')):
            if a==acc:
                try: refs.setdefault(s,[]).append((int(ab),int(b)))
                except: pass
        if not refs: return None
        grp=L('_atom_site.group_PDB'); comp=L('_atom_site.label_comp_id'); ach=L('_atom_site.auth_asym_id'); aseq=L('_atom_site.auth_seq_id'); el=L('_atom_site.type_symbol')
        X=np.array([L('_atom_site.Cartn_x'),L('_atom_site.Cartn_y'),L('_atom_site.Cartn_z')],dtype=float).T
        mod=L('_atom_site.pdbx_PDB_model_num'); first=mod[0]
        keep=np.array([m==first for m in mod]); 
        lig={}
        for i in np.where(keep)[0]:
            if grp[i]=='HETATM' and comp[i] not in EXCL and comp[i]!='HOH' and el[i]!='H': lig.setdefault((comp[i],ach[i],aseq[i]),[]).append(i)
        best=None
        for (c,ch,sq),ix in lig.items():
            if len(ix)<12: continue
            P=[i for i in np.where(keep)[0] if grp[i]=='ATOM' and ach[i] in refs and el[i]!='H']
            if not P: continue
            dist=np.linalg.norm(X[P][:,None,:]-X[ix][None,:,:],axis=2).min(axis=1)
            near=[P[j] for j in np.where(dist<=5.0)[0]]
            if len(set((ach[i],aseq[i]) for i in near))<8: continue
            if best is None or len(near)>best[2]: best=(c,near,len(near))
        if best is None: return None
        res=set()
        for i in best[1]:
            for ab,b in refs[ach[i]]:
                try: res.add(int(aseq[i])-ab+b)
                except: pass
        return {'lig':best[0],'res':sorted(res)}
    return cached(f,f'pocket_{pdb}_{acc}.json',f)
def target_pocket(acc):
    def f():
        ids=search(acc,2.8,['X-RAY DIFFRACTION'])
        if len(ids)<1: ids=search(acc,3.6,['X-RAY DIFFRACTION','ELECTRON MICROSCOPY'])
        out=[];seen=set()
        for p in ids[:25]:
            r=pocket(p,acc)
            if r and r['lig'] not in seen: seen.add(r['lig']); out.append(dict(pdb=p,**r))
            if len(out)>=5: break
        if not out: return {'n':0,'res':[]}
        cnt={}
        for o in out:
            for x in o['res']: cnt[x]=cnt.get(x,0)+1
        thr=2 if len(out)>=3 else 1
        return {'n':len(out),'pdbs':[o['pdb'] for o in out],'res':sorted(x for x,c in cnt.items() if c>=thr)}
    return cached(f,f'tp_{acc}.json',f)
M=substitution_matrices.load('BLOSUM62'); al=PairwiseAligner(); al.substitution_matrix=M; al.open_gap_score=-10; al.extend_gap_score=-0.5; al.mode='global'
def divergence(sa,sb,ra,rb):
    aln=al.align(sa,sb)[0]; ca={};cb={}
    for (a0,a1),(b0,b1) in zip(*aln.aligned):
        for k in range(a1-a0): ca[a0+k+1]=(b0+k+1)
    inv={v:k for k,v in ca.items()}
    UA={ca[r] for r in ra if r in ca}; UB={r for r in rb}
    cols=set()
    for r in ra: cols.add(('a',r))
    # columns: use A-coordinates; B residues mapped to A
    colsA={r for r in ra}|{inv[r] for r in rb if r in inv}
    ident=sum(1 for c in colsA if c in ca and sa[c-1]==sb[ca[c]-1])
    seqid=sum(1 for c in ca if sa[c-1]==sb[ca[c]-1])/max(len(sa),len(sb))
    return 1-ident/max(len(colsA),1),1-seqid
if __name__=='__main__':
    pairs=[l.strip().split('/') for l in open('data/pairs.txt')]
    genes=sorted({g for p in pairs for g in p})
    with ThreadPoolExecutor(8) as ex: acc=dict(zip(genes,ex.map(accession,genes)))
    accs=sorted({a for a in acc.values() if a})
    with ThreadPoolExecutor(8) as ex: list(ex.map(seq,accs))
    with ThreadPoolExecutor(4) as ex: tp=dict(zip(accs,ex.map(target_pocket,accs)))
    json.dump({'acc':acc,'targets':{a:{k:v for k,v in tp[a].items() if k!='res'}|{'n_res':len(tp[a]['res'])} for a in tp}},open(OUT+'/pocket_targets.json','w'),indent=1)
    rows=[]
    for A,B in pairs:
        a,b=acc.get(A),acc.get(B)
        if not a or not b or not tp.get(a,{}).get('res') or not tp.get(b,{}).get('res'): continue
        sa,sb=seq(a),seq(b)
        if not sa or not sb: continue
        pd_,sd_=divergence(sa,sb,tp[a]['res'],tp[b]['res'])
        rows.append(dict(pair=f'{A}/{B}',pocket_div=pd_,seq_div=sd_,n_struct_A=tp[a]['n'],n_struct_B=tp[b]['n']))
    pd.DataFrame(rows).to_csv(OUT+'/pocket_divergence_by_pair.csv',index=False); print(len(rows),'pairs with pocket data')
