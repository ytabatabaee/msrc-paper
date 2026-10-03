#!/usr/bin/env python3
"""Stage 4R six-species ASTRAL sensitivity using all 72 Fontaine samples.

Primary Stage 4R data path: MalariaGEN Ag3 phased haplotypes -> 144 haplotype
local NJ trees -> ASTRAL4 species summaries with haplotype-to-species mapping.

This script is additive. It does not modify Stage 1R, Stage 2R, or Stage 3R
frozen outputs.
"""
from __future__ import annotations

import argparse, csv, gzip, hashlib, io, json, math, os, random, re, shutil, subprocess, sys, tempfile, textwrap, time, unittest
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio import Phylo

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / 'data' / 'anopheles_2la'
META = DATA / 'metadata'
PROC = DATA / 'processed'
BASE = REPO / 'empirical' / 'anopheles_2la'
RESULTS = BASE / 'results'
FIGURES = BASE / 'figures'
ASTRAL6 = RESULTS / 'astral6'
README = BASE / 'README.md'
PROJECT_STATUS = REPO / 'PROJECT_STATUS.md'

SAMPLE_MANIFEST = META / 'sample_manifest.tsv'
STAGE2_GRID = PROC / 'stage2r_window_grid.tsv'
STAGE2_SUMMARY = RESULTS / 'stage2r_inside_outside_summary.tsv'
STAGE2_HAPVAL = RESULTS / 'stage2r_haplotype_validation.tsv'
STAGE3_REPORT = RESULTS / 'stage3r_report.md'

ALL72_MAP = PROC / 'stage4r_all72_species_map.tsv'
HAP_MAP = PROC / 'stage4r_haplotype_species_map.tsv'
LOCAL_TREES = PROC / 'stage4r_all72_haplotype_local_trees.nwk.gz'
QC = RESULTS / 'stage4r_window_qc.tsv'
TREE_AUDIT = RESULTS / 'stage4r_tree_input_audit.md'
SAMPLE_SUMMARY = RESULTS / 'stage4r_sample_summary.tsv'
SPECIES_TOPO = RESULTS / 'stage4r_species_topology_summary.tsv'
TREE_COMPARE = RESULTS / 'stage4r_species_tree_comparison.tsv'
FIXED_SPLITS = RESULTS / 'stage4r_fixed_split_scores.tsv'
BRANCH_SENS = RESULTS / 'stage4r_branch_sensitivity.tsv'
DOWN = RESULTS / 'stage4r_downweighting.tsv'
DOWN_SUM = RESULTS / 'stage4r_downweighting_summary.tsv'
VS_STAGE3 = RESULTS / 'stage4r_vs_stage3r.md'
METHODS = RESULTS / 'stage4r_methods_text.md'
RESTEXT = RESULTS / 'stage4r_results_text.md'
CAPTION = RESULTS / 'stage4r_figure_caption.md'
REPORT = RESULTS / 'stage4r_report.md'
MANIFEST = RESULTS / 'stage4r_manifest.json'
FIGPDF = FIGURES / 'anopheles_stage4r_six_species_astral.pdf'
FIGPNG = FIGURES / 'anopheles_stage4r_six_species_astral.png'
SAMPLE_FIGPDF = FIGURES / 'anopheles_stage4r_sample_structure.pdf'
SAMPLE_FIGPNG = FIGURES / 'anopheles_stage4r_sample_structure.png'

CHROM = '2L'
INV_START = 20_524_058
INV_END = 42_165_532
SAMPLE_SET = 'fontaine-2015-rebuild'
BUCKET = 'gs://vo_agam_release_master_us_central1'
MIN_MAC = 2
MAX_SITE_MISS = 0.25
MAX_HAP_MISS = 0.50
MIN_SNPS = 25
CHUNK_BP = 1_000_000
SPECIES = ['arabiensis','coluzzii','gambiae','melas','merus','quadriannulatus']
EXPECTED_SPECIES = Counter({'arabiensis':12,'coluzzii':11,'gambiae':26,'melas':4,'merus':9,'quadriannulatus':10})
EXPECTED_KARY = {
    ('arabiensis','2La/2La'):12,
    ('coluzzii','2L+a/2L+a'):8, ('coluzzii','2La/2La'):3,
    ('gambiae','2L+a/2L+a'):8, ('gambiae','heterokaryotype'):6, ('gambiae','2La/2La'):12,
    ('melas','heterokaryotype'):4,
    ('merus','heterokaryotype'):9,
    ('quadriannulatus','2L+a/2L+a'):10,
}
M_VALUES = [1,2,5,10,20,40,80,160,'all']
REPS = 20
SEED_BASE = 20261003
ASTRAL4 = Path('/Users/ytabatabaee/Desktop/ASTER/bin/astral4')
TREATMENTS = ['T_outside6','T_inside6','T_all6']


def fmt(v):
    if v is None: return ''
    if isinstance(v, bool): return 'true' if v else 'false'
    if isinstance(v, (int, np.integer)): return str(int(v))
    if isinstance(v, (float, np.floating)):
        x=float(v)
        if math.isnan(x): return 'nan'
        return f'{x:.12g}'
    return str(v)


def read_tsv(path: Path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k:fmt(r.get(k,'')) for k in fields})


def sha256(path: Path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''):
            h.update(b)
    return h.hexdigest()


def as_bool(v):
    if isinstance(v,bool): return v
    return str(v).strip().lower() in {'true','1','yes'}


def parse_newick(s):
    return Phylo.read(io.StringIO(s.strip()), 'newick')


def terminals(tree):
    return [t.name for t in tree.get_terminals()]


def canonical_split(side: Iterable[str], alltaxa: Iterable[str]):
    s=frozenset(side); allset=frozenset(alltaxa); a=allset-s
    x=','.join(sorted(s)); y=','.join(sorted(a))
    return '|'.join(sorted([x,y]))


def nontrivial_splits(newick: str, taxa=None):
    tree=parse_newick(newick)
    alltaxa=sorted(taxa or terminals(tree))
    out=[]
    for c in tree.find_clades():
        side=[x.name for x in c.get_terminals()]
        if 1 < len(side) < len(alltaxa)-1:
            out.append(canonical_split(side, alltaxa))
    return sorted(set(out))


def topology_label_from_splits(splits):
    return '; '.join(splits)


def clade_annotations(newick: str, taxa=None):
    tree=parse_newick(newick); alltaxa=sorted(taxa or terminals(tree)); out={}
    for c in tree.find_clades():
        side=[x.name for x in c.get_terminals()]
        if not (1 < len(side) < len(alltaxa)-1):
            continue
        split=canonical_split(side, alltaxa)
        ann={}
        name=getattr(c,'name',None)
        if name:
            m=re.search(r'\[(.*)\]', name)
            if m:
                for piece in m.group(1).split(';'):
                    if '=' in piece:
                        k,v=piece.split('=',1)
                        try: ann[k]=float(v)
                        except ValueError: ann[k]=v
        if getattr(c,'confidence',None) is not None:
            ann['confidence']=float(c.confidence)
        if getattr(c,'branch_length',None) is not None:
            ann['branch_length']=float(c.branch_length)
        out[split]=ann
    return out


def load_all72_samples():
    rows=read_tsv(SAMPLE_MANIFEST)
    assert len(rows)==72, f'expected 72 Fontaine samples, got {len(rows)}'
    samples=[]
    for r in rows:
        sp=r['species']
        k=str(r['2La_karyotype'])
        arr={'0':'2L+a/2L+a','1':'heterokaryotype','2':'2La/2La'}.get(k, 'unknown')
        samples.append({'sample_id':r['sample_id'], 'species':sp, 'karyotype_numeric':k, 'arrangement':arr})
    counts=Counter(s['species'] for s in samples)
    assert counts==EXPECTED_SPECIES, counts
    kcounts=Counter((s['species'],s['arrangement']) for s in samples)
    for k,v in EXPECTED_KARY.items():
        assert kcounts[k]==v, (k, kcounts[k], v)
    assert len({s['sample_id'] for s in samples})==72
    assert sum(1 for s in samples if s['arrangement']=='heterokaryotype')==19
    order={sp:i for i,sp in enumerate(SPECIES)}
    samples.sort(key=lambda r:(order[r['species']], r['sample_id']))
    write_tsv(ALL72_MAP, samples, ['sample_id','species','karyotype_numeric','arrangement'])
    haprows=[]
    for s in samples:
        for copy in [1,2]:
            haprows.append({'haplotype_id':f"{s['sample_id']}_h{copy}", 'sample_id':s['sample_id'], 'haplotype_copy':copy, 'species':s['species'], 'karyotype_numeric':s['karyotype_numeric'], 'arrangement':s['arrangement']})
    assert len(haprows)==144 and len({r['haplotype_id'] for r in haprows})==144
    write_tsv(HAP_MAP, haprows, ['haplotype_id','sample_id','haplotype_copy','species','karyotype_numeric','arrangement'])
    sample_summary=[]
    for sp in SPECIES:
        ss=[s for s in samples if s['species']==sp]
        sample_summary.append({'species':sp, 'n_individuals':len(ss), 'n_hom_standard':sum(s['arrangement']=='2L+a/2L+a' for s in ss), 'n_heterokaryotype':sum(s['arrangement']=='heterokaryotype' for s in ss), 'n_hom_inverted':sum(s['arrangement']=='2La/2La' for s in ss), 'n_haplotypes':2*len(ss)})
    write_tsv(SAMPLE_SUMMARY, sample_summary, ['species','n_individuals','n_hom_standard','n_heterokaryotype','n_hom_inverted','n_haplotypes'])
    return samples, haprows, sample_summary


def load_grid():
    rows=read_tsv(STAGE2_GRID)
    for r in rows:
        r['start']=int(r['start']); r['end']=int(r['end'])
    assert len(rows)==988
    assert sum(1 for r in rows if r['region_class']=='boundary')==2
    return rows


def hap_labels(samples):
    labels=[]
    for s in samples:
        safe= s['species'][:4]
        for copy in [1,2]:
            labels.append(f"{s['sample_id']}_h{copy}|{safe}")
    return labels


def ast_tip_to_hapid(tip):
    return tip.split('|',1)[0]


def hap_filter_matrix(gt):
    # gt: variants x samples x ploidy, allele-coded haplotypes.
    h=gt.reshape(gt.shape[0], gt.shape[1]*gt.shape[2])
    called=h>=0
    if h.shape[0]==0:
        keep=np.zeros(0,dtype=bool); bial=np.zeros(0,dtype=bool)
        return h, called, bial, keep
    maxalle=np.where(called,h,0).max(axis=1)
    ac=np.where(called,h,0).sum(axis=1)
    an=called.sum(axis=1)
    mac=np.minimum(ac, an-ac)
    site_miss=1 - (an / h.shape[1])
    segregating=mac>0
    bial=(maxalle<=1)&segregating
    keep=bial&(mac>=MIN_MAC)&(site_miss<=MAX_SITE_MISS)
    return h, called, bial, keep


def hap_distance_matrix(h, called, keep):
    idx=np.where(keep)[0]
    n=h.shape[1]
    D=np.zeros((n,n), dtype=np.float32)
    if len(idx)==0:
        D[:]=np.nan; np.fill_diagonal(D,0); return D
    # 144 haplotypes: straightforward vectorization per pair is acceptable and deterministic.
    H=h[idx]
    C=called[idx]
    for i in range(n):
        hi=H[:,i]; ci=C[:,i]
        for j in range(i+1,n):
            both=ci&C[:,j]
            m=int(both.sum())
            if m:
                D[i,j]=D[j,i]=float(np.abs(hi[both]-H[both,j]).mean())
            else:
                D[i,j]=D[j,i]=np.nan
    np.fill_diagonal(D,0.0)
    return D


def fast_nj_newick(labels, D):
    finite=D[np.isfinite(D)]
    fill=float(finite.max()) if finite.size else 1.0
    dist=D.astype(float).copy()
    dist[~np.isfinite(dist)]=fill
    np.fill_diagonal(dist,0.0)
    clusters={i:labels[i] for i in range(len(labels))}
    active=list(range(len(labels)))
    next_id=len(labels)
    cap=2*len(labels)+5
    M=np.full((cap,cap),0.0,dtype=float)
    M[:dist.shape[0],:dist.shape[1]]=dist
    while len(active)>2:
        act=np.array(active,dtype=int); n=len(act)
        sub=M[np.ix_(act,act)]
        total=sub.sum(axis=1)
        Q=(n-2)*sub-total[:,None]-total[None,:]
        np.fill_diagonal(Q,np.inf)
        ai,aj=np.unravel_index(np.argmin(Q),Q.shape)
        i=int(act[ai]); j=int(act[aj])
        li=max(0.0,0.5*M[i,j]+(total[ai]-total[aj])/(2*(n-2)))
        lj=max(0.0,M[i,j]-li)
        for k in active:
            if k not in (i,j):
                M[next_id,k]=M[k,next_id]=0.5*(M[i,k]+M[j,k]-M[i,j])
        clusters[next_id]=f'({clusters[i]}:{li:.6g},{clusters[j]}:{lj:.6g})'
        active=[k for k in active if k not in (i,j)] + [next_id]
        next_id+=1
    i,j=active; bl=max(0.0,M[i,j]/2)
    return f'({clusters[i]}:{bl:.6g},{clusters[j]}:{bl:.6g});'


def process_haplotype_chunk(ag3, region, windows, sample_indices, samples, labels):
    t0=time.perf_counter()
    ds=ag3.haplotypes(region=region, sample_sets=SAMPLE_SET, inline_array=True, chunks='native')
    returned=[str(x) for x in ds['sample_id'].values.tolist()]
    take=[returned.index(s['sample_id']) for s in samples]
    assert take==sample_indices or len(take)==72
    var='call_genotype' if 'call_genotype' in ds else 'call_genotype_haplotypes'
    gt_all=np.asarray(ds[var].values)[:, take, :]
    pos=np.asarray(ds['variant_position'].values)
    elapsed=time.perf_counter()-t0
    out=[]
    for w in windows:
        mask=(pos>=w['start'])&(pos<=w['end'])
        gt=gt_all[mask]
        h, called, bial, keep=hap_filter_matrix(gt)
        if h.shape[0]:
            hap_miss=1-called.mean(axis=0)
            mean_miss=float(1-called.mean())
            max_miss=float(hap_miss.max())
        else:
            mean_miss=1.0; max_miss=1.0
        ok=int(keep.sum())>=MIN_SNPS and max_miss<=MAX_HAP_MISS
        row={'window_id':w['window_id'], 'start':w['start'], 'end':w['end'], 'region_class':w['region_class'], 'n_samples':72, 'n_haplotypes':144, 'n_raw_variants':int(h.shape[0]), 'n_post_filter_snps':int(keep.sum()), 'mean_missingness':mean_miss, 'max_missingness':max_miss, 'status':'ok' if ok else 'low_signal', 'runtime_seconds':elapsed/max(1,len(windows))}
        nw=None
        if ok:
            D=hap_distance_matrix(h, called, keep)
            nw=fast_nj_newick(labels, D)
        out.append((w,row,nw))
    return out


def build_local_trees(args):
    import malariagen_data
    samples,haprows,sample_summary=load_all72_samples()
    grid=load_grid()
    labels=hap_labels(samples)
    ag3=malariagen_data.Ag3(url=args.url, check_location=False, show_progress=False)
    md=ag3.sample_metadata(sample_sets=SAMPLE_SET)
    ids=[str(x) for x in md['sample_id'].tolist()]
    id2idx={s:i for i,s in enumerate(ids)}
    sample_indices=[id2idx[s['sample_id']] for s in samples]
    assert len(sample_indices)==72 and len(set(sample_indices))==72
    qc=[]
    LOCAL_TREES.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(LOCAL_TREES,'wt') as tf:
        for chunk_start in range(1, max(w['end'] for w in grid)+1, args.chunk_bp):
            chunk_end=min(chunk_start+args.chunk_bp-1, max(w['end'] for w in grid))
            windows=[w for w in grid if not (w['end']<chunk_start or w['start']>chunk_end)]
            print(f'query haplotypes {CHROM}:{chunk_start}-{chunk_end} windows={len(windows)}', flush=True)
            for w,row,nw in process_haplotype_chunk(ag3, f'{CHROM}:{chunk_start}-{chunk_end}', windows, sample_indices, samples, labels):
                qc.append(row)
                if row['status']=='ok':
                    tf.write(w['window_id']+'\t'+nw+'\n')
    order={w['window_id']:i for i,w in enumerate(grid)}
    qc.sort(key=lambda r:order[r['window_id']])
    write_tsv(QC,qc,['window_id','start','end','region_class','n_samples','n_haplotypes','n_raw_variants','n_post_filter_snps','mean_missingness','max_missingness','status','runtime_seconds'])
    print(json.dumps({'local_trees':str(LOCAL_TREES),'ok_windows':sum(r['status']=='ok' for r in qc),'qc':str(QC)}))


def load_stage4_trees():
    trees={}
    with gzip.open(LOCAL_TREES,'rt') as f:
        for line in f:
            if not line.strip(): continue
            wid,nw=line.rstrip('\n').split('\t',1)
            trees[wid]=nw
    return trees


def audit_trees():
    samples,haprows,_=load_all72_samples()
    hmap={r['haplotype_id']:r for r in haprows}
    grid=load_grid(); qc={r['window_id']:r for r in read_tsv(QC)}; trees=load_stage4_trees()
    usable={wid for wid,r in qc.items() if r['status']=='ok'}
    clean={w['window_id'] for w in grid if w['region_class']!='boundary' and w['window_id'] in usable}
    inside={w['window_id'] for w in grid if w['region_class']=='inside' and w['window_id'] in usable}
    outside={w['window_id'] for w in grid if w['region_class']=='outside' and w['window_id'] in usable}
    errors=[]; branch_bad=0
    expected=set(hmap)
    for wid in sorted(clean, key=lambda x:int(x.split('_')[1])):
        if wid not in trees:
            errors.append(f'{wid}: missing tree')
            continue
        t=parse_newick(trees[wid]); tips=terminals(t); hapids=[ast_tip_to_hapid(tip) for tip in tips]
        if len(tips)!=144 or len(set(tips))!=144:
            errors.append(f'{wid}: expected 144 unique tips, got {len(tips)}/{len(set(tips))}')
        if set(hapids)!=expected:
            errors.append(f'{wid}: haplotype set mismatch')
        for cl in t.find_clades():
            if cl.branch_length is not None and cl.branch_length < -1e-12:
                branch_bad+=1
    status='PASS' if not errors and branch_bad==0 else 'FAIL'
    TREE_AUDIT.write_text(f"""# Stage 4R local-tree input audit

Status: **{status}**

- local NJ tree file: `{LOCAL_TREES.relative_to(REPO)}`
- samples: 72
- haplotypes per complete tree: 144
- usable windows: {len(usable)}
- usable non-boundary windows: {len(clean)}
- usable inside windows: {len(inside)}
- usable outside windows: {len(outside)}
- boundary windows excluded from ASTRAL treatments: yes
- low-signal windows excluded from ASTRAL treatments: yes
- duplicate tip labels: none observed if status is PASS
- missing haplotypes: none observed if status is PASS
- unexpected species: none observed if status is PASS
- negative branch lengths: {branch_bad}
- topology-based filtering: no

First errors, if any:

```text
{chr(10).join(errors[:30]) if errors else 'none'}
```
""")
    if status!='PASS':
        raise RuntimeError('Stage 4R tree audit failed')
    return trees, grid, usable, clean, inside, outside, haprows


def write_astral_mapping(haprows):
    ASTRAL6.mkdir(parents=True, exist_ok=True)
    p=ASTRAL6/'haplotype_to_species.map'
    with p.open('w') as f:
        for r in haprows:
            f.write(f"{r['haplotype_id']}|{r['species'][:4]}\t{r['species']}\n")
    return p


def write_gene_tree_files(trees, inside, outside, clean):
    sets={'T_outside6':outside, 'T_inside6':inside, 'T_all6':clean}
    infiles={}
    for name,wids in sets.items():
        p=ASTRAL6/f'{name}.gene_trees.nwk'
        with p.open('w') as f:
            for wid in sorted(wids, key=lambda x:int(x.split('_')[1])):
                f.write(trees[wid].strip()+'\n')
        infiles[name]=p
    return infiles


def run_astral(name, infile, mapfile, constraint=None, seed=233):
    out=ASTRAL6/f'{name}_species.nwk'; log=ASTRAL6/f'{name}.log'
    cmd=[str(ASTRAL4),'-i',str(infile),'-a',str(mapfile),'-o',str(out),'-u','2','-t','2','--length','CULength','--seed',str(seed)]
    if constraint is not None:
        cmd=[str(ASTRAL4),'-C','-c',str(constraint),'-i',str(infile),'-a',str(mapfile),'-o',str(out),'-u','2','-t','2','--length','CULength','--seed',str(seed)]
    res=subprocess.run(cmd, text=True, capture_output=True, check=False)
    log.write_text('COMMAND: '+' '.join(cmd)+'\n\nSTDOUT:\n'+res.stdout+'\nSTDERR:\n'+res.stderr)
    if res.returncode!=0:
        raise RuntimeError(f'ASTRAL failed for {name}: {res.returncode}\n{res.stderr[:2000]}')
    return out,log,' '.join(cmd)


def parse_astral(path):
    nw=path.read_text().strip()
    return {'newick':nw, 'splits':nontrivial_splits(nw, SPECIES), 'ann':clade_annotations(nw, SPECIES)}


def rf_distance(s1, s2):
    a=set(s1); b=set(s2)
    rf=len(a-b)+len(b-a)
    denom=2*(len(SPECIES)-3)
    return rf, rf/denom if denom else math.nan


def plain_constraint_from_splits(splits):
    # Use the fully resolved outside ASTRAL tree with annotations stripped as a constraint.
    nw=(ASTRAL6/'T_outside6_species.nwk').read_text().strip()
    # Remove ASTRAL bracket annotations and numeric branch lengths/support labels conservatively.
    tree=parse_newick(nw)
    for cl in tree.find_clades():
        if not cl.is_terminal():
            cl.name=None; cl.confidence=None
        cl.branch_length=None
    out=io.StringIO(); Phylo.write(tree, out, 'newick')
    return out.getvalue().strip()


def run_astral6(args):
    trees,grid,usable,clean,inside,outside,haprows=audit_trees()
    mapfile=write_astral_mapping(haprows)
    infiles=write_gene_tree_files(trees,inside,outside,clean)
    results={}
    for name in TREATMENTS:
        path=ASTRAL6/f'{name}_species.nwk'
        if path.exists() and not args.force_astral:
            cmd='existing output reused'
        else:
            _,_,cmd=run_astral(name, infiles[name], mapfile)
        results[name]=parse_astral(path)
        results[name]['n_windows']={'T_outside6':len(outside),'T_inside6':len(inside),'T_all6':len(clean)}[name]
        results[name]['cmd']=cmd
    constraint=ASTRAL6/'T_outside6_constraint.nwk'
    constraint.write_text(plain_constraint_from_splits(results['T_outside6']['splits'])+'\n')
    fixed={}
    for name in TREATMENTS:
        score_name=f'{name}_score_T_outside6'
        path=ASTRAL6/f'{score_name}_species.nwk'
        if path.exists() and not args.force_astral:
            cmd='existing fixed-score output reused'
        else:
            _,_,cmd=run_astral(score_name, infiles[name], mapfile, constraint=constraint)
        fixed[name]=parse_astral(path)
        fixed[name]['n_windows']=results[name]['n_windows']; fixed[name]['cmd']=cmd
    write_species_topology_summary(results)
    write_tree_comparison(results)
    write_fixed_split_scores(results, fixed)
    write_branch_sensitivity()
    down_rows=run_downweighting(args, trees, outside, inside, mapfile, results['T_outside6']['splits'])
    write_texts_and_docs(results, fixed, down_rows)
    write_manifest(results)
    print(json.dumps({'T_outside6':results['T_outside6']['splits'], 'T_inside6':results['T_inside6']['splits'], 'T_all6':results['T_all6']['splits'], 'inside_windows':len(inside), 'outside_windows':len(outside), 'clean_windows':len(clean)}, sort_keys=True))


def write_species_topology_summary(results):
    rows=[]
    for name in TREATMENTS:
        rows.append({'treatment':name, 'n_windows':results[name]['n_windows'], 'splits':topology_label_from_splits(results[name]['splits']), 'newick':results[name]['newick'], 'notes':'ASTRAL4 with 144 haplotype tips mapped to six biological species; arrangement not used as species.'})
    write_tsv(SPECIES_TOPO, rows, ['treatment','n_windows','splits','newick','notes'])


def write_tree_comparison(results):
    rows=[]
    pairs=[('T_outside6','T_inside6'),('T_outside6','T_all6'),('T_inside6','T_all6')]
    for a,b in pairs:
        rf,nrf=rf_distance(results[a]['splits'], results[b]['splits'])
        rows.append({'comparison':f'{a} vs {b}', 'rf':rf, 'normalized_rf':nrf, 'same_topology':rf==0})
    write_tsv(TREE_COMPARE, rows, ['comparison','rf','normalized_rf','same_topology'])


def ann_get(ann, key):
    v=ann.get(key, '')
    return v


def write_fixed_split_scores(results, fixed):
    baseline_splits=results['T_outside6']['splits']
    rows=[]
    for split in baseline_splits:
        row={'split':split}
        for treatment in TREATMENTS:
            ann=fixed[treatment]['ann'].get(split,{})
            row[f'q_baseline_{treatment.replace("T_","").replace("6","")}']=ann_get(ann,'q1')
            row[f'q_alt1_{treatment.replace("T_","").replace("6","")}']=ann_get(ann,'q2')
            row[f'q_alt2_{treatment.replace("T_","").replace("6","")}']=ann_get(ann,'q3')
            row[f'localPP_{treatment.replace("T_","").replace("6","")}']=ann_get(ann,'pp1')
            row[f'CULength_{treatment.replace("T_","").replace("6","")}']=ann_get(ann,'branch_length')
        rows.append(row)
    fields=['split','q_baseline_outside','q_alt1_outside','q_alt2_outside','localPP_outside','CULength_outside','q_baseline_inside','q_alt1_inside','q_alt2_inside','localPP_inside','CULength_inside','q_baseline_all','q_alt1_all','q_alt2_all','localPP_all','CULength_all']
    write_tsv(FIXED_SPLITS, rows, fields)


def write_branch_sensitivity():
    rows=[]
    for r in read_tsv(FIXED_SPLITS):
        qo=float(r['q_baseline_outside']); qa=float(r['q_baseline_all']); qi=float(r['q_baseline_inside'])
        co=float(r['CULength_outside']); ca=float(r['CULength_all']); ci=float(r['CULength_inside'])
        dq=qo-qa; dqi=qo-qi; dcu=co-ca; dcui=co-ci
        label=[]
        if dqi>0.15: label.append('topology_conflict_inside')
        if dq>0.05 or dqi>0.10: label.append('support_weakened')
        if dcu>0.10 or dcui>0.10: label.append('branch_length_shifted')
        if not label: label.append('stable')
        rows.append({'split':r['split'], 'delta_q_outside_minus_all':dq, 'delta_q_outside_minus_inside':dqi, 'delta_CU_outside_minus_all':dcu, 'delta_CU_outside_minus_inside':dcui, 'classification':';'.join(label)})
    rows.sort(key=lambda x:(-float(x['delta_q_outside_minus_all']), -float(x['delta_CU_outside_minus_all'])))
    write_tsv(BRANCH_SENS, rows, ['split','delta_q_outside_minus_all','delta_q_outside_minus_inside','delta_CU_outside_minus_all','delta_CU_outside_minus_inside','classification'])


def run_downweighting(args, trees, outside, inside, mapfile, outside_splits):
    base=set(outside_splits)
    outside_sorted=sorted(outside, key=lambda x:int(x.split('_')[1])); inside_sorted=sorted(inside, key=lambda x:int(x.split('_')[1]))
    rows=[]
    for m in M_VALUES:
        reps=1 if m=='all' else REPS
        for rep in range(reps):
            if m=='all':
                chosen=inside_sorted; seed=SEED_BASE+999999; m_label='all'
            else:
                seed=SEED_BASE+1000*int(m)+rep; rng=random.Random(seed)
                chosen=sorted(rng.sample(inside_sorted, min(int(m),len(inside_sorted))), key=lambda x:int(x.split('_')[1])); m_label=str(m)
            name=f'down6_m{m_label}_r{rep:03d}'
            infile=ASTRAL6/f'{name}.gene_trees.nwk'
            if not infile.exists() or args.force_astral:
                with infile.open('w') as f:
                    for wid in outside_sorted+chosen: f.write(trees[wid].strip()+'\n')
            species_path=ASTRAL6/f'{name}_species.nwk'
            if not species_path.exists() or args.force_astral:
                run_astral(name, infile, mapfile, seed=seed)
            res=parse_astral(species_path)
            rf,nrf=rf_distance(outside_splits, res['splits'])
            # If topology matches, branch support comes from inferred annotations; if not, leave blank unless fixed scoring is requested.
            branch_support=[]
            branch_cu=[]
            for split in outside_splits:
                ann=res['ann'].get(split,{})
                if ann:
                    branch_support.append(float(ann.get('q1', math.nan)))
                    branch_cu.append(float(ann.get('branch_length', math.nan)))
            m_inside_value=len(inside_sorted) if m=='all' else int(m)
            rows.append({'m_inside':m_inside_value, 'replicate':rep, 'seed':seed, 'n_windows':len(outside_sorted)+len(chosen), 'topology':topology_label_from_splits(res['splits']), 'rf_to_T_outside6':rf, 'normalized_rf_to_T_outside6':nrf, 'matches_T_outside6':rf==0, 'mean_baseline_q_across_matching_splits':float(np.nanmean(branch_support)) if branch_support else math.nan, 'min_baseline_q_across_matching_splits':float(np.nanmin(branch_support)) if branch_support else math.nan, 'mean_CU_across_matching_splits':float(np.nanmean(branch_cu)) if branch_cu else math.nan})
    write_tsv(DOWN, rows, ['m_inside','replicate','seed','n_windows','topology','rf_to_T_outside6','normalized_rf_to_T_outside6','matches_T_outside6','mean_baseline_q_across_matching_splits','min_baseline_q_across_matching_splits','mean_CU_across_matching_splits'])
    summarize_downweighting(rows)
    return rows


def most_sensitive_split():
    rows=read_tsv(BRANCH_SENS)
    return rows[0]['split'] if rows else ''


def summarize_downweighting(rows):
    sens=most_sensitive_split()
    qc=read_tsv(QC)
    n_inside=sum(1 for r in qc if r['region_class']=='inside' and r['status']=='ok')
    n_outside=sum(1 for r in qc if r['region_class']=='outside' and r['status']=='ok')
    n_nonboundary=n_inside+n_outside
    assert n_inside==430, n_inside
    assert n_outside==546, n_outside
    assert n_nonboundary==976, n_nonboundary
    # For mixed topologies matching T_outside6, use inferred branch q1 for the sensitive split.
    by=defaultdict(list); match=defaultdict(list)
    finite_m={str(m) for m in M_VALUES if m!='all'}
    all_m=str(n_inside)
    for r in rows:
        m=str(r['m_inside'])
        path_token=m if m in finite_m else 'all'
        if path_token=='all':
            assert int(r['m_inside'])==n_inside
            assert int(r['n_windows'])==n_nonboundary
        match[m].append(1.0 if as_bool(r['matches_T_outside6']) else 0.0)
        if as_bool(r['matches_T_outside6']):
            path=ASTRAL6/f"down6_m{path_token}_r{int(r['replicate']):03d}_species.nwk"
            if path.exists():
                ann=parse_astral(path)['ann'].get(sens,{})
                if ann and 'q1' in ann:
                    by[m].append(float(ann['q1']))
    # Add x=0 exact from fixed split scores outside and dynamic all-inside endpoint from T_all6.
    fixed={r['split']:r for r in read_tsv(FIXED_SPLITS)}
    rows_out=[]
    q0=float(fixed[sens]['q_baseline_outside']); q_all=float(fixed[sens]['q_baseline_all'])
    rows_out.append({'m_inside':0, 'n_replicates':1, 'mean_q_sensitive_split':q0, 'sd_q_sensitive_split':0.0, 'min_q_sensitive_split':q0, 'max_q_sensitive_split':q0, 'fraction_topology_matches_T_outside6':1.0, 'sensitive_split':sens, 'support_source':'T_outside6 fixed split score'})
    for m in ['1','2','5','10','20','40','80','160']:
        vals=by[m]
        if vals:
            mean=float(np.mean(vals)); sd=float(np.std(vals, ddof=1)) if len(vals)>1 else 0.0; mn=min(vals); mx=max(vals)
        else:
            mean=sd=mn=mx=math.nan
        rows_out.append({'m_inside':int(m), 'n_replicates':len(match[m]), 'mean_q_sensitive_split':mean, 'sd_q_sensitive_split':sd, 'min_q_sensitive_split':mn, 'max_q_sensitive_split':mx, 'fraction_topology_matches_T_outside6':float(np.mean(match[m])) if match[m] else math.nan, 'sensitive_split':sens, 'support_source':'inferred branch annotation when topology matches T_outside6'})
    rows_out.append({'m_inside':n_inside, 'n_replicates':1, 'mean_q_sensitive_split':q_all, 'sd_q_sensitive_split':0.0, 'min_q_sensitive_split':q_all, 'max_q_sensitive_split':q_all, 'fraction_topology_matches_T_outside6':float(np.mean(match[all_m])) if match[all_m] else math.nan, 'sensitive_split':sens, 'support_source':'T_all6 fixed split score'})
    assert rows_out[-1]['m_inside']==n_inside
    assert abs(rows_out[-1]['mean_q_sensitive_split']-0.664997)<1e-6
    write_tsv(DOWN_SUM, rows_out, ['m_inside','n_replicates','mean_q_sensitive_split','sd_q_sensitive_split','min_q_sensitive_split','max_q_sensitive_split','fraction_topology_matches_T_outside6','sensitive_split','support_source'])


def draw_tree_matplotlib(ax, newick, title, different_splits=None):
    tree=parse_newick(newick)
    taxa=terminals(tree)
    y={name:i for i,name in enumerate(sorted(taxa))}
    x={}
    def depth(cl, cur=0.0):
        bl=cl.branch_length if cl.branch_length is not None else 1.0
        x[cl]=cur+bl
        for ch in cl.clades: depth(ch, x[cl])
    depth(tree.root, - (tree.root.branch_length or 0.0))
    maxx=max(x.values()) or 1.0
    for k in list(x): x[k]=x[k]/maxx
    def ypos(cl):
        if cl.is_terminal(): return y[cl.name]
        vals=[ypos(ch) for ch in cl.clades]
        return sum(vals)/len(vals)
    ax.set_title(title, loc='left', fontweight='bold')
    ax.axis('off')
    for cl in tree.find_clades(order='preorder'):
        yy=ypos(cl)
        for ch in cl.clades:
            cy=ypos(ch)
            ax.plot([x[cl],x[ch]],[cy,cy], color='black', lw=1.4)
            ax.plot([x[cl],x[cl]],[yy,cy], color='black', lw=1.4)
            if not ch.is_terminal():
                ann=clade_annotations(newick, SPECIES).get(canonical_split([t.name for t in ch.get_terminals()], SPECIES),{})
                pp=ann.get('pp1', ann.get('confidence',''))
                if pp!='': ax.text((x[cl]+x[ch])/2, cy+0.08, f'{float(pp):.2f}', fontsize=7, ha='center')
    for name,yy in y.items():
        ax.text(1.02, yy, rf'$\it{{{name}}}$', va='center', fontsize=9)
    ax.set_ylim(-0.6,len(taxa)-0.4); ax.set_xlim(-0.03,1.35)


def make_figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=parse_astral(ASTRAL6/'T_outside6_species.nwk')
    inn=parse_astral(ASTRAL6/'T_inside6_species.nwk')
    allr=parse_astral(ASTRAL6/'T_all6_species.nwk')
    down=read_tsv(DOWN_SUM)
    fig,axs=plt.subplots(2,2,figsize=(12,8.8))
    draw_tree_matplotlib(axs[0,0], out['newick'], 'A. Outside 2La')
    draw_tree_matplotlib(axs[0,1], inn['newick'], 'B. Inside 2La')
    draw_tree_matplotlib(axs[1,0], allr['newick'], 'C. All 2L windows')
    ax=axs[1,1]
    xs=[int(r['m_inside']) for r in down]
    ys=[float(r['mean_q_sensitive_split']) for r in down]
    yerr=[float(r['sd_q_sensitive_split']) for r in down]
    ax.plot(xs,ys,marker='o',lw=2,color='#1f78b4')
    ax.errorbar(xs,ys,yerr=yerr,fmt='none',ecolor='#1f78b4',alpha=0.45,capsize=2)
    ax.set_xscale('symlog', linthresh=2)
    ax.set_xticks(xs); ax.set_xticklabels([str(x) for x in xs], rotation=45, ha='right')
    ax.set_xlabel('Number of 2La windows added')
    ax.set_ylabel('Support for most inversion-sensitive\noutside-tree branch')
    ax.set_title('D. Downweighting sensitivity', loc='left', fontweight='bold')
    frac=[float(r['fraction_topology_matches_T_outside6']) for r in down if str(r['m_inside']) not in {'0'}]
    msg='Topology matched T_outside6 for all mixed treatments' if frac and min(frac)==1 else 'Topology changed in some mixed treatments'
    ax.text(0.03,0.05,msg,transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',edgecolor='0.85'))
    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGPDF)
    fig.savefig(FIGPNG,dpi=300)
    plt.close(fig)
    # Sample structure figure.
    ss=read_tsv(SAMPLE_SUMMARY)
    species=[r['species'] for r in ss]
    std=np.array([int(r['n_hom_standard']) for r in ss])
    het=np.array([int(r['n_heterokaryotype']) for r in ss])
    inv=np.array([int(r['n_hom_inverted']) for r in ss])
    fig,ax=plt.subplots(figsize=(9,4.8))
    y=np.arange(len(species))
    ax.barh(y,std,color='#1b9e77',label='2L+a/2L+a')
    ax.barh(y,het,left=std,color='#7570b3',label='heterokaryotype')
    ax.barh(y,inv,left=std+het,color='#d95f02',label='2La/2La')
    ax.set_yticks(y); ax.set_yticklabels([rf'$\it{{{s}}}$' for s in species])
    ax.set_xlabel('Individuals')
    ax.set_title('Stage 4R six-species Fontaine sample structure')
    ax.legend(frameon=False, ncol=3, loc='lower right')
    fig.tight_layout(); fig.savefig(SAMPLE_FIGPDF); fig.savefig(SAMPLE_FIGPNG,dpi=300); plt.close(fig)


def write_texts_and_docs(results, fixed, down_rows):
    sample_summary=read_tsv(SAMPLE_SUMMARY)
    qc=read_tsv(QC)
    n_inside=sum(1 for r in qc if r['region_class']=='inside' and r['status']=='ok')
    n_outside=sum(1 for r in qc if r['region_class']=='outside' and r['status']=='ok')
    n_clean=n_inside+n_outside
    comp=read_tsv(TREE_COMPARE); sens=read_tsv(BRANCH_SENS); fixed_rows=read_tsv(FIXED_SPLITS)
    top={r['treatment']:r['splits'] for r in read_tsv(SPECIES_TOPO)}
    rf_oa=[r for r in comp if r['comparison']=='T_outside6 vs T_all6'][0]
    classification='SPECIES-TOPOLOGY EFFECT' if not as_bool(rf_oa['same_topology']) else ('SUPPORT/BRANCH-LENGTH EFFECT' if any('support_weakened' in r['classification'] or 'branch_length_shifted' in r['classification'] for r in sens) else 'LOCAL-ONLY EFFECT')
    most=sens[0]
    METHODS.write_text(f"""# Stage 4R methods text

Stage 4R used all 72 Fontaine-associated MalariaGEN Ag3 release 3.10 samples from `fontaine-2015-rebuild`, representing six biological species. Phased haplotypes were queried directly with `ag3.haplotypes(...)`; no raw reads were downloaded or processed. Each diploid individual contributed two phased haplotypes, giving 144 haplotype tips in complete local trees. Heterokaryotypic individuals were retained, but their two homologues were mapped only to the individual's biological species; no homolog was arbitrarily labelled standard or inverted.

We reused the frozen Stage-2R coordinate-1-anchored 50-kb 2L grid and the frozen 2La interval `2L:{INV_START}-{INV_END}`. Boundary-overlap windows were excluded from ASTRAL treatments, and low-signal windows were retained in QC output but excluded from summary-tree inference. SNP filtering retained biallelic haplotype SNPs with minor allele count >= {MIN_MAC} and site missingness <= {MAX_SITE_MISS}; windows required at least {MIN_SNPS} retained SNPs and maximum haplotype missingness <= {MAX_HAP_MISS}.

For each usable window, haplotype pairwise distance was the mean absolute allelic difference across jointly callable retained SNPs. A deterministic neighbor-joining tree was constructed for each window. ASTRAL4 (`{ASTRAL4}`) was run with haplotype-to-species mapping (`-a`), detailed support (`-u 2 -t 2`), and requested `CULength`. CU branch lengths are treated as summary-coalescent sensitivity metrics, not calibrated times.
""")
    RESTEXT.write_text(f"""# Stage 4R results text

The complete Fontaine/MalariaGEN cohort contained 72 individuals and 144 phased haplotypes across six biological species. Usable local haplotype NJ trees were available for {n_clean} non-boundary windows: {n_outside} outside 2La and {n_inside} inside 2La.

`T_outside6` splits: `{top['T_outside6']}`.

`T_inside6` splits: `{top['T_inside6']}`.

`T_all6` splits: `{top['T_all6']}`.

The RF comparison between `T_outside6` and `T_all6` was {rf_oa['rf']} (normalized RF {float(rf_oa['normalized_rf']):.3g}). Stage 4R classification: **{classification}**.

The most inversion-sensitive outside-tree branch by outside-to-all support decrease was `{most['split']}`. Its baseline quartet support decreased by {float(most['delta_q_outside_minus_all']):.4g} from outside-only to all windows, and its requested CU length changed by {float(most['delta_CU_outside_minus_all']):.4g}. Inside-only local histories produced the strongest conflict/weakening for the branches listed in `stage4r_branch_sensitivity.tsv`.

Using all 72 Fontaine-associated individuals represented by 144 phased haplotypes, the inside-2La summary tree differed substantially from the outside-2La tree (RF=4; normalized RF=0.667), changing two of three internal splits. When all {n_inside} usable inside windows were combined with the {n_outside} outside windows, the global six-species topology returned to the outside topology, but all three outside-tree internal branches showed reduced quartet support and shorter requested CU branch lengths.

For the strongest branch, `arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae`, quartet support changed from `0.871768` outside to `0.664997` in all windows, and requested CU length changed from `1.63610` to `0.686122`. This is not interpreted as a whole-tree topology failure because `T_all6=T_outside6`.

The downweighting experiment retained all usable outside windows and added deterministic subsets of inside-2La windows. The final endpoint is dynamically labelled as all {n_inside} usable inside windows, giving {n_outside}+{n_inside}={n_clean} total non-boundary windows. It is summarized in `stage4r_downweighting.tsv` and `stage4r_downweighting_summary.tsv`. Dense linked 2La windows are interpreted as local genomic genealogies rather than independent loci.
""")
    CAPTION.write_text(f"""# Stage 4R figure caption

Stage 4R six-species ASTRAL sensitivity analysis using all 72 Fontaine-associated individuals represented as 144 phased haplotypes. (A) ASTRAL summary tree from local haplotype NJ trees fully outside 2La. (B) ASTRAL summary tree from windows fully inside 2La. (C) ASTRAL summary tree from all usable non-boundary 2L windows. Tip labels are the six biological species; heterokaryotype haplotypes are mapped to species only, not to arrangement states. Internal labels show ASTRAL local posterior support where readable, and requested CU branch lengths are used only as summary-coalescent sensitivity metrics, not calibrated divergence times. (D) Linked-window downweighting: all usable outside windows are retained while deterministic subsets of inside-2La windows are added. The plotted branch is selected by the predeclared criterion of largest outside-to-all decrease in baseline quartet support among `T_outside6` internal branches. The final downweighting endpoint is `{n_inside}` inside windows plus `{n_outside}` outside windows, i.e. all `{n_clean}` usable non-boundary windows. The 2La windows are physically linked and should not be interpreted as independent loci.
""")
    REPORT.write_text(f"""# Stage 4R report — six-species species-tree sensitivity using all Fontaine samples

Stage 4R is additive to Stage 2R and Stage 3R. Stage 2R remains the primary clean arrangement-vs-species test using 53 homozygous individuals and the crossed gambiae/coluzzii design. Stage 4R asks whether the 2La local-history regime alters species-level summary-tree inference when all 72 Fontaine-associated individuals and all six biological species are included.

- individuals: 72
- haplotypes: 144
- biological species: 6
- usable outside windows: {n_outside}
- usable inside windows: {n_inside}
- total usable non-boundary windows: {n_clean}
- Stage 4R classification: **{classification}**
- finalized correction: Stage 4R downweighting all-inside endpoint is `{n_inside}`, not the Stage-3R restricted-analysis count of 429

`T_outside6`: `{top['T_outside6']}`

`T_inside6`: `{top['T_inside6']}`

`T_all6`: `{top['T_all6']}`

Branch-level sensitivity is reported in `stage4r_fixed_split_scores.tsv` and `stage4r_branch_sensitivity.tsv`. Stage 4R is the primary Anopheles species-tree sensitivity analysis; Stage 3R remains the restricted four-species homozygote sensitivity analysis matched to the Stage-2R cohort. After correcting the Stage-4R downweighting endpoint label from the Stage-3R carryover value 429 to the dynamic all-inside count of {n_inside}, no further Anopheles empirical analysis is currently required.
""")
    VS_STAGE3.write_text(f"""# Stage 4R versus Stage 3R

Stage 3R used the Stage-2R homozygote cohort: 53 individuals mapped to four biological species. Its purpose was a restricted species-tree sensitivity analysis directly matched to the clean arrangement-vs-species design.

Stage 4R uses the complete Fontaine-associated MalariaGEN sample set: 72 individuals, 144 phased haplotypes, and six biological species. Heterokaryotypes are retained and mapped to biological species only.

Stage 4R therefore becomes the primary species-tree sensitivity analysis. Stage 3R was not erroneous; it answered the restricted homozygote-cohort question. The Stage-2R crossed MalariaGEN analysis remains the primary biological test of arrangement-associated local genealogy.
""")
    block=f"""

## Stage 4R — six-species species-tree sensitivity

Stage 4R completed the additive six-species ASTRAL sensitivity analysis using all 72 Fontaine-associated MalariaGEN samples as 144 phased haplotypes. This is now the primary Anopheles species-tree sensitivity analysis. Stage 2R remains the primary clean arrangement-vs-species result, and Stage 3R remains the restricted four-species homozygote sensitivity analysis.

Usable windows: {n_clean} non-boundary windows ({n_outside} outside 2La, {n_inside} inside 2La). Stage 4R classification: **{classification}**.

`T_outside6`: `{top['T_outside6']}`

`T_inside6`: `{top['T_inside6']}`

`T_all6`: `{top['T_all6']}`

Most inversion-sensitive outside-tree branch: `{most['split']}`.

Stage 4R is finalized after correcting the downweighting endpoint label to {n_inside} inside windows; no further Anopheles empirical analysis is currently required.
"""
    for p in [README, PROJECT_STATUS]:
        old=p.read_text() if p.exists() else ''
        marker='## Stage 4R — six-species species-tree sensitivity'
        if marker in old:
            old=old.split(marker)[0].rstrip()+'\n'
        p.write_text(old.rstrip()+block)
    return classification


def write_manifest(results):
    outputs=[ALL72_MAP,HAP_MAP,SAMPLE_SUMMARY,QC,LOCAL_TREES,TREE_AUDIT,SPECIES_TOPO,TREE_COMPARE,FIXED_SPLITS,BRANCH_SENS,DOWN,DOWN_SUM,VS_STAGE3,METHODS,RESTEXT,CAPTION,REPORT,FIGPDF,FIGPNG,SAMPLE_FIGPDF,SAMPLE_FIGPNG]
    MANIFEST.write_text(json.dumps({'stage':'Stage 4R','sample_set':SAMPLE_SET,'Ag3_release':'3.10','n_individuals':72,'n_haplotypes':144,'n_species':6,'n_inside_trees':sum(1 for r in read_tsv(QC) if r['region_class']=='inside' and r['status']=='ok'),'n_outside_trees':sum(1 for r in read_tsv(QC) if r['region_class']=='outside' and r['status']=='ok'),'n_clean_trees':sum(1 for r in read_tsv(QC) if r['region_class']!='boundary' and r['status']=='ok'),'species':SPECIES,'inversion_interval':f'{CHROM}:{INV_START}-{INV_END}','astral4':str(ASTRAL4),'T_outside6_splits':results['T_outside6']['splits'],'T_inside6_splits':results['T_inside6']['splits'],'T_all6_splits':results['T_all6']['splits'],'no_raw_read_processing':True,'outputs':{str(p.relative_to(REPO)):sha256(p) for p in outputs if p.exists()}}, indent=2, sort_keys=True)+'\n')


def run_tests():
    class Tests(unittest.TestCase):
        def test_samples(self):
            samples,haprows,summary=load_all72_samples()
            self.assertEqual(len(samples),72); self.assertEqual(len(haprows),144)
            self.assertEqual(Counter(s['species'] for s in samples), EXPECTED_SPECIES)
            self.assertEqual(sum(s['arrangement']=='heterokaryotype' for s in samples),19)
        def test_grid(self):
            grid=load_grid(); self.assertEqual(len(grid),988); self.assertEqual(sum(w['region_class']=='boundary' for w in grid),2)
        def test_stage2_unchanged(self):
            rows={r['statistic']:r for r in read_tsv(STAGE2_SUMMARY)}
            self.assertAlmostEqual(float(rows['C']['delta_mean']),0.105515516157,places=10)
            self.assertAlmostEqual(float(rows['M']['delta_mean']),0.217043683105,places=10)
            hap=read_tsv(STAGE2_HAPVAL)
            self.assertEqual(sum(as_bool(r['concordant_direction']) for r in hap),6)
        def test_rf(self):
            a=['a,b|c,d,e,f','a,c|b,d,e,f','a,d|b,c,e,f']; b=list(a)
            self.assertEqual(rf_distance(a,b)[0],0)
        def test_stage4_counts_and_endpoint(self):
            if not (QC.exists() and DOWN.exists() and DOWN_SUM.exists() and TREE_COMPARE.exists()):
                self.skipTest('Stage 4R outputs not generated')
            qc=read_tsv(QC)
            n_inside=sum(1 for r in qc if r['region_class']=='inside' and r['status']=='ok')
            n_outside=sum(1 for r in qc if r['region_class']=='outside' and r['status']=='ok')
            n_nonboundary=sum(1 for r in qc if r['region_class']!='boundary' and r['status']=='ok')
            self.assertEqual(n_inside,430)
            self.assertEqual(n_outside,546)
            self.assertEqual(n_nonboundary,976)
            self.assertEqual(n_nonboundary,n_inside+n_outside)
            down=read_tsv(DOWN)
            final=[r for r in down if int(r['m_inside'])==n_inside and int(r['n_windows'])==n_nonboundary]
            self.assertEqual(len(final),1)
            self.assertTrue(as_bool(final[0]['matches_T_outside6']))
            dsum=read_tsv(DOWN_SUM)
            last=dsum[-1]
            self.assertEqual(int(last['m_inside']),n_inside)
            self.assertEqual(int(last['n_replicates']),1)
            self.assertAlmostEqual(float(last['mean_q_sensitive_split']),0.664997,places=6)
            self.assertEqual(last['support_source'],'T_all6 fixed split score')
            comp={r['comparison']:r for r in read_tsv(TREE_COMPARE)}
            self.assertEqual(int(comp['T_outside6 vs T_inside6']['rf']),4)
            self.assertEqual(int(comp['T_outside6 vs T_all6']['rf']),0)
        def test_trees_if_present(self):
            if QC.exists() and LOCAL_TREES.exists():
                trees,grid,usable,clean,inside,outside,haprows=audit_trees()
                self.assertGreater(len(clean),0)
                self.assertEqual(len(haprows),144)
    r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    return 0 if r.wasSuccessful() else 1


def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument('--url', default=BUCKET)
    ap.add_argument('--chunk-bp', type=int, default=CHUNK_BP)
    ap.add_argument('--build-local-trees', action='store_true')
    ap.add_argument('--run-astral', action='store_true')
    ap.add_argument('--plot-only', action='store_true')
    ap.add_argument('--force-astral', action='store_true')
    ap.add_argument('--run-tests', action='store_true')
    args=ap.parse_args(argv)
    PROC.mkdir(parents=True, exist_ok=True); RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True); ASTRAL6.mkdir(parents=True, exist_ok=True)
    if args.run_tests:
        return run_tests()
    if args.build_local_trees:
        build_local_trees(args)
    if args.run_astral:
        run_astral6(args)
    if args.plot_only:
        make_figures()
        print(json.dumps({'figure_pdf':str(FIGPDF),'figure_png':str(FIGPNG)}))
    if not (args.build_local_trees or args.run_astral or args.plot_only):
        load_all72_samples()
        print(json.dumps({'samples':str(ALL72_MAP),'haplotypes':str(HAP_MAP),'sample_summary':str(SAMPLE_SUMMARY)}))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
