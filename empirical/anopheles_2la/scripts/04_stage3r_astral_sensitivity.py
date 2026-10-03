#!/usr/bin/env python3
"""Stage 3R ASTRAL/ASTER species-tree sensitivity and Stage-2R reporting fixes."""
from __future__ import annotations

import argparse, csv, gzip, hashlib, io, json, math, os, random, re, shutil, subprocess, sys, tempfile, textwrap, unittest
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'msrc-paper-mpl'))
os.environ.setdefault('XDG_CACHE_HOME', str(Path(tempfile.gettempdir()) / 'msrc-paper-xdg'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo

REPO=Path(__file__).resolve().parents[3]
DATA=REPO/'data'/'anopheles_2la'
PROC=DATA/'processed'
BASE=REPO/'empirical'/'anopheles_2la'
RESULTS=BASE/'results'
FIGURES=BASE/'figures'
ASTRAL_DIR=RESULTS/'astral'
README=BASE/'README.md'
PROJECT_STATUS=REPO/'PROJECT_STATUS.md'

PRIMARY=PROC/'stage1r_primary_homozygote_samples.tsv'
GRID=PROC/'stage2r_window_grid.tsv'
QC=RESULTS/'stage2r_window_qc.tsv'
TREES=PROC/'stage2r_local_nj_trees.nwk.gz'
CROSS=RESULTS/'stage2r_crossed_distance_signal.tsv'
MARGIN=RESULTS/'stage2r_quartet_margin_signal.tsv'
TOPO=RESULTS/'stage2r_quartet_topology_signal.tsv'
SUMMARY=RESULTS/'stage2r_inside_outside_summary.tsv'
CIRC=RESULTS/'stage2r_circular_shift_tests.tsv'
PHYS=RESULTS/'stage2r_physical_interval_tests.tsv'
HAPVAL=RESULTS/'stage2r_haplotype_validation.tsv'
STAGE2_RESULTS=RESULTS/'stage2r_results_text.md'
STAGE2_CAPTION=RESULTS/'stage2r_main_figure_caption.md'
STAGE2_REPORT=RESULTS/'stage2r_report.md'
STAGE2_MANIFEST=RESULTS/'stage2r_manifest.json'
MAIN_PDF=FIGURES/'anopheles_2la_main.pdf'
MAIN_PNG=FIGURES/'anopheles_2la_main.png'
STAGE2_FULL_PDF=FIGURES/'stage2r_full_2L_spatial_signal.pdf'
STAGE2_FULL_PNG=FIGURES/'stage2r_full_2L_spatial_signal.png'

MAP=PROC/'stage3r_individual_species_map.tsv'
AUDIT=RESULTS/'stage3r_tree_input_audit.md'
ENV=RESULTS/'stage3r_astral_environment.md'
TOPO_SUM=RESULTS/'stage3r_species_topology_summary.tsv'
FIXED=RESULTS/'stage3r_fixed_topology_scores.tsv'
DOWN=RESULTS/'stage3r_downweighting.tsv'
DOWN_SUMMARY=RESULTS/'stage3r_downweighting_summary.tsv'
BRANCH=RESULTS/'stage3r_branch_length_summary.tsv'
ARRCTRL=RESULTS/'stage3r_arrangement_class_quartet_summary.tsv'
METHODS=RESULTS/'stage3r_methods_text.md'
RESTEXT=RESULTS/'stage3r_results_text.md'
CAPTION=RESULTS/'stage3r_figure_caption.md'
REPORT=RESULTS/'stage3r_report.md'
MANIFEST=RESULTS/'stage3r_manifest.json'
FIG=FIGURES/'anopheles_stage3r_astral_sensitivity.pdf'
FIGPNG=FIGURES/'anopheles_stage3r_astral_sensitivity.png'

INV_START_MB=20.524058; INV_END_MB=42.165532
SPECIES=['arabiensis','coluzzii','gambiae','quadriannulatus']
T_OUTSIDE_SPLIT='arabiensis,quadriannulatus|coluzzii,gambiae'
T_INSIDE_SPLIT='arabiensis,gambiae|coluzzii,quadriannulatus'
T_ALT_SPLIT='arabiensis,coluzzii|gambiae,quadriannulatus'
ASTRAL4=Path('/Users/ytabatabaee/Desktop/ASTER/bin/astral4')
M_VALUES=[1,2,5,10,20,40,80,160,'all']
SEED_BASE=20261003
REPS=20  # ASTRAL runs are external processes; 20 deterministic reps keeps this sensitivity cheap and reproducible.
EXPECTED_STAGE2={
    ('C','mean_inside'):0.0960668024032, ('C','mean_outside'):-0.00944871375335, ('C','delta_mean'):0.105515516157,
    ('M','mean_inside'):0.196433856604, ('M','mean_outside'):-0.0206098265003, ('M','delta_mean'):0.217043683105,
}

def read_tsv(p:Path):
    with p.open(newline='') as f: return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(p:Path, rows, fields):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', newline='') as f:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator='\n'); w.writeheader()
        for r in rows: w.writerow({k:fmt(r.get(k,'')) for k in fields})

def fmt(v):
    if v is None: return ''
    if isinstance(v,bool): return 'true' if v else 'false'
    if isinstance(v,float): return 'nan' if math.isnan(v) else f'{v:.12g}'
    return str(v)

def as_bool(v):
    if isinstance(v,bool): return v
    return str(v).strip().lower() in {'true','1','yes'}

def sha256(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def parse_newick(s): return Phylo.read(io.StringIO(s.strip()), 'newick')

def terminals(tree): return [t.name for t in tree.get_terminals()]

def canonical_split(side:Iterable[str], alltaxa:Iterable[str]):
    s=frozenset(side); a=frozenset(alltaxa)-s
    x=','.join(sorted(s)); y=','.join(sorted(a))
    return '|'.join(sorted([x,y]))

def tree_topology(newick:str):
    t=parse_newick(newick); taxa=sorted(terminals(t)); splits=[]
    for c in t.find_clades():
        side=[x.name for x in c.get_terminals()]
        if 1<len(side)<len(taxa)-1: splits.append(canonical_split(side,taxa))
    if not splits and len(taxa)==4:
        # Bio may collapse support label strangely; find clade with two terminals.
        for c in t.find_clades():
            side=[x.name for x in c.get_terminals()]
            if len(side)==2: splits.append(canonical_split(side,taxa))
    return sorted(set(splits))[0] if splits else 'unresolved'

def topology_label(split):
    if split=='unresolved': return split
    a,b=split.split('|')
    return f'({a})|({b})'

def clade_annotations(newick):
    out={}; t=parse_newick(newick); taxa=sorted(terminals(t))
    for c in t.find_clades():
        side=[x.name for x in c.get_terminals()]
        if not (1<len(side)<len(taxa)-1): continue
        split=canonical_split(side,taxa)
        ann={}
        name=getattr(c,'name',None)
        if name:
            m=re.search(r'\[(.*)\]', name)
            if m:
                for piece in m.group(1).split(';'):
                    if '=' in piece:
                        k,v=piece.split('=',1)
                        try: ann[k]=float(v)
                        except ValueError: pass
        if getattr(c,'confidence',None) is not None: ann['confidence']=float(c.confidence)
        if getattr(c,'branch_length',None) is not None: ann['branch_length']=float(c.branch_length)
        out[split]=ann
    return out

def load_stage2():
    grid=read_tsv(GRID); qc=read_tsv(QC)
    g={r['window_id']:r for r in grid}; q={r['window_id']:r for r in qc}
    usable={wid for wid,r in q.items() if r['status']=='ok'}
    clean={wid for wid in usable if g[wid]['region_class']!='boundary'}
    inside={wid for wid in clean if g[wid]['region_class']=='inside'}; outside={wid for wid in clean if g[wid]['region_class']=='outside'}
    return grid,qc,inside,outside,clean

def fix_stage2_texts():
    hap=read_tsv(HAPVAL); n=sum(1 for r in hap if as_bool(r['concordant_direction'])); total=len(hap)
    assert n==6 and total==6
    for p in [STAGE2_RESULTS, README, PROJECT_STATUS, STAGE2_REPORT]:
        if not p.exists(): continue
        s=p.read_text()
        s=s.replace('concordant in direction for ' + (chr(48) + '/6') + ' windows','concordant in direction for 6/6 windows')
        s=s.replace('Haplotype validation was concordant for ' + (chr(48) + '/6') + ' deterministic validation windows.','Haplotype validation was concordant for 6/6 deterministic validation windows.')
        s=s.replace('Haplotype validation was concordant for ' + (chr(48) + '/6') + ' deterministic validation windows','Haplotype validation was concordant for 6/6 deterministic validation windows')
        p.write_text(s)
    # Patch generator robustly for future runs.
    script=BASE/'scripts'/'03_stage2r_full_scan.py'
    if script.exists():
        s=script.read_text()
        s=s.replace("sum(1 for r in hap_rows if r['concordant_direction']=='true')", "sum(1 for r in hap_rows if str(r['concordant_direction']).lower() == 'true' or r['concordant_direction'] is True)")
        script.write_text(s)

def check_stage2_numbers():
    rows={r['statistic']:r for r in read_tsv(SUMMARY)}
    for (stat,key),exp in EXPECTED_STAGE2.items():
        got=float(rows[stat][key]); assert abs(got-exp)<1e-10, (stat,key,got,exp)

def make_main_figure():
    c=read_tsv(CROSS); m=read_tsv(MARGIN); d=read_tsv(TOPO); q=read_tsv(QC)
    maps={'C':{r['window_id']:r for r in c}, 'M':{r['window_id']:r for r in m}, 'D':{r['window_id']:r for r in d}, 'Q':{r['window_id']:r for r in q}}
    grid=read_tsv(GRID); summary={r['statistic']:r for r in read_tsv(SUMMARY)}; circ={r['statistic']:r for r in read_tsv(CIRC)}
    fig, axes=plt.subplots(4,1,figsize=(10,8.8),sharex=True,gridspec_kw={'height_ratios':[1,1,1,0.75]})
    panels=[('C','C','C(w)'),('M','mean_M','M(w)'),('D','D',r'D(w) = q$_{arrangement}$ - q$_{species}$'),('Q','n_post_filter_snps','Retained SNPs')]
    for ax,(name,key,ylabel) in zip(axes,panels):
        xs=[]; ys=[]; colors=[]; markers=[]
        for w in grid:
            wid=w['window_id']; val=maps[name][wid][key]
            try: y=float(val)
            except Exception: y=math.nan
            xs.append((int(w['start'])+int(w['end']))/2/1e6); ys.append(y)
            rc=w['region_class']; colors.append('#d95f02' if rc=='inside' else '#757575' if rc=='boundary' else '#1b9e77'); markers.append('s' if rc=='boundary' else 'o')
        ax.axvspan(INV_START_MB, INV_END_MB, color='#756bb1', alpha=0.16, lw=0)
        ax.axhline(0,color='0.45',lw=0.8,ls='--')
        for mk in ['o','s']:
            sel=[i for i,z in enumerate(markers) if z==mk and not math.isnan(ys[i])]
            ax.scatter([xs[i] for i in sel],[ys[i] for i in sel],s=7 if mk=='o' else 16,c=[colors[i] for i in sel],marker=mk,linewidths=0,alpha=0.9)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.tick_params(labelsize=9)
        if name in {'C','M'}:
            txt=f"inside mean = {float(summary[name]['mean_inside']):.3g}\noutside mean = {float(summary[name]['mean_outside']):.3g}\ncircular p = {float(circ[name]['p_one_sided_ge_observed']):.3g}"
            ax.text(0.012,0.96,txt, transform=ax.transAxes, va='top', ha='left', fontsize=8.5, bbox=dict(facecolor='white', edgecolor='0.85', alpha=0.9, pad=2.5))
        ax.text((INV_START_MB+INV_END_MB)/2, ax.get_ylim()[1], '2La', ha='center', va='top', fontsize=9, color='#4b3f8f')
    axes[-1].set_xlabel('Genomic position on 2L (Mb)', fontsize=11)
    axes[0].set_title('Anopheles 2La arrangement-vs-species signal across chromosome arm 2L', fontsize=12)
    fig.tight_layout()
    fig.savefig(MAIN_PDF)
    fig.savefig(MAIN_PNG, dpi=300)
    # Also replace full spatial signal with same axis-based figure for consistency.
    fig.savefig(STAGE2_FULL_PDF)
    fig.savefig(STAGE2_FULL_PNG, dpi=300)
    plt.close(fig)
    STAGE2_CAPTION.write_text(f"""# Stage 2R main figure caption

Full chromosome-arm 2L scan of the 2La arrangement-associated genealogy signal. Points are unsmoothed 50-kb windows plotted by genomic midpoint in Mb. The purple shaded region marks the frozen 2La interval, 20.524058-42.165532 Mb. Boundary-overlap windows are shown with square markers and excluded from primary inside-vs-outside tests. Panels show: (A) crossed gambiae/coluzzii distance contrast `C(w)`, (B) continuous crossed-quartet margin `M(w)`, (C) discrete quartet topology support `D(w)=q_arrangement-q_species`, and (D) retained post-filter SNP count. Horizontal dashed lines mark zero for the signal panels. C and M annotations report inside and outside means and exact circular-shift p-values.
""")

def load_trees():
    trees={}
    with gzip.open(TREES,'rt') as f:
        for line in f:
            if not line.strip(): continue
            wid,nw=line.rstrip('\n').split('\t',1); trees[wid]=nw
    return trees

def write_species_map():
    rows=[r for r in read_tsv(PRIMARY) if r['include_primary']=='true']
    rows.sort(key=lambda r:r['sample_id'])
    out=[{'sample_id':r['sample_id'],'species':r['species'],'arrangement':r['arrangement']} for r in rows]
    counts=Counter(r['species'] for r in out); assert counts==Counter({'arabiensis':12,'coluzzii':11,'gambiae':20,'quadriannulatus':10})
    assert not any(r['arrangement']=='heterokaryotype' for r in out)
    write_tsv(MAP,out,['sample_id','species','arrangement'])
    return out

def tip_to_sample(tip): return tip.split('|',1)[0]

def audit_trees(trees, clean, inside, outside, mapping):
    map_ids={r['sample_id'] for r in mapping}; errors=[]; branch_bad=0
    for wid in clean:
        nw=trees.get(wid)
        if not nw: errors.append(f'missing tree {wid}'); continue
        t=parse_newick(nw); tips=terminals(t)
        if len(tips)!=53 or len(set(tips))!=53: errors.append(f'{wid}: expected 53 unique tips got {len(tips)}/{len(set(tips))}')
        if {tip_to_sample(x) for x in tips}!=map_ids: errors.append(f'{wid}: tip sample set mismatch')
        for cl in t.find_clades():
            if cl.branch_length is not None and cl.branch_length < -1e-12: branch_bad+=1
    status='PASS' if not errors and branch_bad==0 else 'FAIL'
    AUDIT.write_text(f"""# Stage 3R tree input audit

Status: **{status}**

- local NJ tree file: `{TREES.relative_to(REPO)}`
- usable non-boundary trees: {len(clean)}
- inside usable trees: {len(inside)}
- outside usable trees: {len(outside)}
- expected tips per tree: 53
- duplicate-tip errors: none observed if status is PASS
- missing-individual errors: none observed if status is PASS
- negative branch lengths: {branch_bad}
- boundary windows excluded: yes
- low-signal windows excluded: yes

First errors, if any:

```text
{chr(10).join(errors[:20]) if errors else 'none'}
```
""")
    if status!='PASS': raise RuntimeError('tree audit failed')

def astral_help():
    out=subprocess.run([str(ASTRAL4),'-h'],text=True,capture_output=True,check=False)
    return out.stdout+out.stderr

def write_env(help_text):
    sha=sha256(ASTRAL4) if ASTRAL4.exists() else 'missing'
    commit=subprocess.run(['git','-C',str(ASTRAL4.parents[1]),'rev-parse','HEAD'],text=True,capture_output=True).stdout.strip() if ASTRAL4.exists() else ''
    ENV.write_text(f"""# Stage 3R ASTRAL/ASTER environment

- executable: `{ASTRAL4}`
- executable present: `{ASTRAL4.exists()}`
- executable SHA256: `{sha}`
- version/help first line: `{help_text.splitlines()[0] if help_text.splitlines() else 'unavailable'}`
- ASTER git commit: `{commit}`
- direct multi-individual mapping support: yes, help documents `-a/--mapping`, a gene-name to taxon-name mapping file with one gene name followed by one taxon name separated by space or tab.
- primary command template: `astral4 -i <trees> -a <mapping> -o <output> -u 2 -t 2 --length CULength`
- fixed scoring template: `astral4 -C -c <outside_tree> -i <trees> -a <mapping> -o <output> -u 2 -t 2 --length CULength`

## Help excerpt

```text
{help_text[:4000]}
```
""")

def prepare_astral_inputs(trees, inside, outside, clean, mapping):
    ASTRAL_DIR.mkdir(parents=True, exist_ok=True)
    mapfile=ASTRAL_DIR/'individual_to_species.map'
    with mapfile.open('w') as f:
        for r in mapping:
            # local tree tips include sample|species-prefix|arrangement tags
            matches=[]
            for wid,nw in list(trees.items())[:1]:
                for tip in terminals(parse_newick(nw)):
                    if tip_to_sample(tip)==r['sample_id']: matches.append(tip)
            if matches: f.write(f'{matches[0]}\t{r["species"]}\n')
    sets={'T_inside':inside,'T_outside':outside,'T_all':clean}
    infiles={}
    for name,wids in sets.items():
        path=ASTRAL_DIR/f'{name}.gene_trees.nwk'
        with path.open('w') as f:
            for wid in sorted(wids, key=lambda x:int(x.split('_')[1])): f.write(trees[wid].strip()+'\n')
        infiles[name]=path
    return mapfile,infiles

def run_astral(name, infile, mapfile, constraint=None, seed=233):
    out=ASTRAL_DIR/f'{name}_species.nwk'; log=ASTRAL_DIR/f'{name}.log'
    cmd=[str(ASTRAL4),'-i',str(infile),'-a',str(mapfile),'-o',str(out),'-u','2','-t','2','--length','CULength','--seed',str(seed)]
    if constraint is not None:
        cmd=[str(ASTRAL4),'-C','-c',str(constraint),'-i',str(infile),'-a',str(mapfile),'-o',str(out),'-u','2','-t','2','--length','CULength','--seed',str(seed)]
    res=subprocess.run(cmd,text=True,capture_output=True,check=False)
    log.write_text('COMMAND: '+' '.join(cmd)+'\n\nSTDOUT:\n'+res.stdout+'\nSTDERR:\n'+res.stderr)
    if res.returncode!=0: raise RuntimeError(f'ASTRAL failed {name}: {res.returncode}\n{res.stderr[:1000]}')
    return out,log,cmd

def parse_support_for_tree(path):
    nw=path.read_text().strip(); split=tree_topology(nw); ann=clade_annotations(nw).get(split,{})
    return nw,split,ann

def all_three_topologies():
    taxa=SPECIES
    return [canonical_split([taxa[0],taxa[1]],taxa), canonical_split([taxa[0],taxa[2]],taxa), canonical_split([taxa[0],taxa[3]],taxa)]



def all_species_splits():
    return {T_OUTSIDE_SPLIT, T_INSIDE_SPLIT, T_ALT_SPLIT}

def canonical_q_map_from_inferred(split, ann):
    """Map ASTRAL q fields to canonical species splits when q1 is the inferred split.

    For q2/q3, ASTRAL's ordering is not guaranteed from the Newick alone, so only
    q1 is assigned unambiguously here. Fixed-topology -C scoring is used whenever
    support for a non-inferred split is required.
    """
    out={k:math.nan for k in all_species_splits()}
    out[split]=float(ann['q1'])
    return out

def fixed_score_rows_by_treatment():
    return {r['treatment']:r for r in read_tsv(FIXED)}

def baseline_quartet_support(treatment=None, split=None, ann=None, fixed_rows=None):
    """Return support specifically for T_OUTSIDE_SPLIT, never raw q1 unless safe.

    If the inferred/scored split is T_OUTSIDE_SPLIT, q1 is the baseline support.
    Otherwise a fixed-topology scoring row must be supplied and q_baseline is used.
    """
    if split == T_OUTSIDE_SPLIT and ann is not None and 'q1' in ann:
        return float(ann['q1']), 'inferred q1 equals T_outside split'
    if treatment is None:
        raise ValueError('treatment is required when inferred split is not T_outside')
    fixed=fixed_rows if fixed_rows is not None else fixed_score_rows_by_treatment()
    row=fixed[treatment]
    return float(row['q_baseline']), f'{treatment} fixed T_outside score'

def load_existing_astral_results():
    paths={'T_outside':ASTRAL_DIR/'T_outside_species.nwk','T_inside':ASTRAL_DIR/'T_inside_species.nwk','T_all':ASTRAL_DIR/'T_all_species.nwk'}
    nwindows={'T_outside':547,'T_inside':429,'T_all':976}
    out={}
    for treatment,path in paths.items():
        nw,split,ann=parse_support_for_tree(path)
        out[treatment]={'tree':nw,'split':split,'ann':ann,'n_windows':nwindows[treatment],'cmd':'existing ASTRAL output reused; no rerun'}
    assert out['T_outside']['split']==T_OUTSIDE_SPLIT
    assert out['T_all']['split']==T_OUTSIDE_SPLIT
    assert out['T_inside']['split']==T_INSIDE_SPLIT
    fixed=fixed_score_rows_by_treatment()
    assert abs(baseline_quartet_support('T_outside', out['T_outside']['split'], out['T_outside']['ann'], fixed)[0]-0.931104)<1e-6
    assert abs(baseline_quartet_support('T_all', out['T_all']['split'], out['T_all']['ann'], fixed)[0]-0.701540)<1e-6
    assert abs(baseline_quartet_support('T_inside', out['T_inside']['split'], out['T_inside']['ann'], fixed)[0]-0.408833)<1e-6
    assert abs(float(out['T_inside']['ann']['q1'])-0.456944)<1e-6
    return out

def species_quartet_scores_from_stage2(wids, crows=None):
    # For four biological species, use all-53 distance class files would not identify species quartet alternatives directly.
    # Here we report ASTRAL annotations where possible; fallback frequencies are not needed because ASTRAL is available.
    return {}


def plain_tree_from_split(split):
    left,right=split.split('|')
    a=left.split(','); b=right.split(',')
    return f"(({a[0]},{a[1]}),({b[0]},{b[1]}));\n"

def do_astral(trees,inside,outside,clean,mapping):
    help_text=astral_help(); write_env(help_text)
    mapfile,infiles=prepare_astral_inputs(trees,inside,outside,clean,mapping)
    results={}; cmds={}
    for t in ['T_outside','T_inside','T_all']:
        out,log,cmd=run_astral(t,infiles[t],mapfile); nw,split,ann=parse_support_for_tree(out); results[t]={'tree':nw,'split':split,'ann':ann,'n_windows':len({'T_outside':outside,'T_inside':inside,'T_all':clean}[t]),'cmd':' '.join(cmd)}; cmds[t]=cmd
    baseline=ASTRAL_DIR/'T_outside_plain_constraint.nwk'
    baseline.write_text(plain_tree_from_split(results['T_outside']['split']))
    fixed=[]
    for t in ['T_outside','T_inside','T_all']:
        out,log,cmd=run_astral(f'{t}_score_T_outside', infiles[t], mapfile, constraint=baseline)
        nw,split,ann=parse_support_for_tree(out)
        fixed.append({'treatment':t,'n_windows':results[t]['n_windows'],'q_baseline':ann.get('q1',ann.get('f1','')),'q_alt1':ann.get('q2',ann.get('f2','')),'q_alt2':ann.get('q3',ann.get('f3','')),'localPP_baseline':ann.get('pp1',ann.get('posterior','')),'topology_scored':topology_label(results['T_outside']['split']),'notes':'ASTRAL4 -C scoring of T_outside topology; support keys preserved when emitted'})
    write_tsv(FIXED,fixed,['treatment','n_windows','q_baseline','q_alt1','q_alt2','localPP_baseline','topology_scored','notes'])
    top_rows=[]
    topo_keys=all_three_topologies()
    for t in ['T_outside','T_inside','T_all']:
        ann=results[t]['ann']
        top_rows.append({'treatment':t,'n_windows':results[t]['n_windows'],'topology':topology_label(results[t]['split']),'localPP_or_support':ann.get('pp1',ann.get('confidence',ann.get('posterior',''))),'quartet1':ann.get('q1',ann.get('f1','')),'quartet2':ann.get('q2',ann.get('f2','')),'quartet3':ann.get('q3',ann.get('f3','')),'notes':'ASTRAL4 multi-individual mapping; arrangement labels not used as species'})
    write_tsv(TOPO_SUM,top_rows,['treatment','n_windows','topology','localPP_or_support','quartet1','quartet2','quartet3','notes'])
    branch=[]
    for t in ['T_outside','T_inside','T_all']:
        ann=results[t]['ann']
        branch.append({'treatment':t,'topology':topology_label(results[t]['split']),'branch_length':ann.get('branch_length',''),'length_type':'CULength requested','notes':'CU lengths are summary-coalescent sensitivity metrics only; not calibrated divergence time'})
    write_tsv(BRANCH,branch,['treatment','topology','branch_length','length_type','notes'])
    return results,fixed,mapfile,infiles

def downweight(trees,inside,outside,results,mapfile):
    outside_sorted=sorted(outside,key=lambda x:int(x.split('_')[1])); inside_sorted=sorted(inside,key=lambda x:int(x.split('_')[1]))
    base_split=results['T_outside']['split']; rows=[]
    for m in M_VALUES:
        reps=1 if m=='all' else REPS
        for rep in range(reps):
            if m=='all': chosen=inside_sorted; seed=SEED_BASE+999999
            else:
                seed=SEED_BASE+1000*int(m)+rep; rng=random.Random(seed); chosen=sorted(rng.sample(inside_sorted, min(int(m),len(inside_sorted))), key=lambda x:int(x.split('_')[1]))
            name=f'down_m{m}_r{rep:03d}'
            infile=ASTRAL_DIR/f'{name}.gene_trees.nwk'
            with infile.open('w') as f:
                for wid in outside_sorted+chosen: f.write(trees[wid].strip()+'\n')
            out,log,cmd=run_astral(name,infile,mapfile,seed=seed)
            nw,split,ann=parse_support_for_tree(out)
            rows.append({'m_inside':m,'replicate':rep,'seed':seed,'n_windows':len(outside_sorted)+len(chosen),'topology':topology_label(split),'matches_T_outside':split==base_split,'support':ann.get('pp1',ann.get('confidence','')),'q1':ann.get('q1',ann.get('f1','')),'q2':ann.get('q2',ann.get('f2','')),'q3':ann.get('q3',ann.get('f3','')),'notes':'ASTRAL4 with all outside windows plus deterministic subset of inside windows'})
    write_tsv(DOWN,rows,['m_inside','replicate','seed','n_windows','topology','matches_T_outside','support','q1','q2','q3','notes'])
    return rows

def arrangement_control():
    rows=[]
    for src,name,key in [(TOPO,'D','D'),(MARGIN,'M','mean_M'),(CROSS,'C','C')]:
        data=read_tsv(src)
        for cls in ['inside','outside']:
            vals=[float(r[key]) for r in data if r['region_class']==cls and str(r[key])!='nan']
            rows.append({'statistic':name,'region_class':cls,'n_windows':len(vals),'mean':sum(vals)/len(vals),'median':sorted(vals)[len(vals)//2],'notes':'arrangement-class quartet/distance diagnostic; these classes are not biological species'})
    write_tsv(ARRCTRL,rows,['statistic','region_class','n_windows','mean','median','notes'])
    return rows

def split_sides(split):
    left,right=split.split('|')
    return left.split(','), right.split(',')

def italic_species(name):
    return rf'$\it{{{name}}}$'

def draw_quartet_tree(ax, split, title, ann, n_windows, extra_note=''):
    left,right=split_sides(split)
    ax.set_title(title, loc='left', fontsize=12, fontweight='bold')
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    x_tip_l, x_cherry_l, x_cherry_r, x_tip_r = 0.07, 0.30, 0.70, 0.93
    y_left=[0.75,0.55]; y_right=[0.45,0.25]
    y_mid_l=sum(y_left)/2; y_mid_r=sum(y_right)/2
    lw=2.0
    ax.plot([x_cherry_l,x_cherry_l],[min(y_left),max(y_left)], color='black', lw=lw)
    ax.plot([x_cherry_r,x_cherry_r],[min(y_right),max(y_right)], color='black', lw=lw)
    for y,taxon in zip(y_left,left):
        ax.plot([x_tip_l,x_cherry_l],[y,y], color='black', lw=lw)
        ax.text(x_tip_l-0.02,y,italic_species(taxon),ha='right',va='center',fontsize=10)
    for y,taxon in zip(y_right,right):
        ax.plot([x_cherry_r,x_tip_r],[y,y], color='black', lw=lw)
        ax.text(x_tip_r+0.02,y,italic_species(taxon),ha='left',va='center',fontsize=10)
    ax.plot([x_cherry_l,x_cherry_r],[y_mid_l,y_mid_r], color='black', lw=lw+0.8)
    pp=float(ann.get('pp1', ann.get('localPP', ann.get('confidence', float('nan')))))
    cu=float(ann.get('branch_length', float('nan')))
    q1=float(ann.get('q1', float('nan'))); q2=float(ann.get('q2', float('nan'))); q3=float(ann.get('q3', float('nan')))
    txt=f"localPP = {pp:.4g}\nCU = {cu:.3f}\nq = ({q1:.3f}, {q2:.3f}, {q3:.3f})\nn = {n_windows} windows"
    if extra_note:
        txt += '\n' + extra_note
    ax.text(0.50,0.96,txt,ha='center',va='top',fontsize=8.5,bbox=dict(facecolor='white', edgecolor='0.82', pad=3))

def baseline_support_summary(down_rows, fixed_rows):
    fixed={r['treatment']:r for r in fixed_rows}
    rows=[]
    outside_q=float(fixed['T_outside']['q_baseline'])
    all_q=float(fixed['T_all']['q_baseline'])
    inside_q=float(fixed['T_inside']['q_baseline'])
    assert abs(outside_q-0.931104)<1e-6
    assert abs(all_q-0.701540)<1e-6
    assert abs(inside_q-0.408833)<1e-6
    rows.append({'m_inside':'0','n_replicates':1,'mean_q_outside_topology':outside_q,'sd_q_outside_topology':0.0,'min_q_outside_topology':outside_q,'max_q_outside_topology':outside_q,'fraction_topology_matches_outside':1.0,'support_source':'T_outside exact'})
    by=defaultdict(list); matches=defaultdict(list); sources=defaultdict(set)
    for r in down_rows:
        m=str(r['m_inside'])
        if m == 'all':
            continue
        inferred_split = T_OUTSIDE_SPLIT if as_bool(r['matches_T_outside']) else None
        if inferred_split == T_OUTSIDE_SPLIT:
            q_base=float(r['q1'])
            source='inferred q1 because topology == T_outside'
        else:
            raise RuntimeError(f"Mixed treatment {m} replicate {r['replicate']} does not match T_outside; fixed scoring is required before summarising.")
        by[m].append(q_base)
        matches[m].append(1.0 if as_bool(r['matches_T_outside']) else 0.0)
        sources[m].add(source)
    for m in ['1','2','5','10','20','40','80','160']:
        vals=by[m]
        assert vals, f'missing downweighting values for m={m}'
        mean=sum(vals)/len(vals)
        sd=(sum((v-mean)**2 for v in vals)/(len(vals)-1))**0.5 if len(vals)>1 else 0.0
        frac=sum(matches[m])/len(matches[m])
        assert abs(frac-1.0)<1e-12, f'mixed treatment m={m} did not always match T_outside'
        rows.append({'m_inside':m,'n_replicates':len(vals),'mean_q_outside_topology':mean,'sd_q_outside_topology':sd,'min_q_outside_topology':min(vals),'max_q_outside_topology':max(vals),'fraction_topology_matches_outside':frac,'support_source':'; '.join(sorted(sources[m]))})
    rows.append({'m_inside':'429','n_replicates':1,'mean_q_outside_topology':all_q,'sd_q_outside_topology':0.0,'min_q_outside_topology':all_q,'max_q_outside_topology':all_q,'fraction_topology_matches_outside':1.0,'support_source':'T_all exact'})
    assert abs(rows[-1]['mean_q_outside_topology']-inside_q)>0.1, 'x=429 must use T_all, not inside-only baseline support'
    write_tsv(DOWN_SUMMARY, rows, ['m_inside','n_replicates','mean_q_outside_topology','sd_q_outside_topology','min_q_outside_topology','max_q_outside_topology','fraction_topology_matches_outside','support_source'])
    return rows

def make_stage3_fig(results,down_rows,fixed_rows=None):
    fixed_rows = fixed_rows if fixed_rows is not None else read_tsv(FIXED)
    summary_rows=baseline_support_summary(down_rows, fixed_rows)
    fig,axs=plt.subplots(2,2,figsize=(11,8.2))
    draw_quartet_tree(axs[0,0], results['T_outside']['split'], 'A. Outside 2La', results['T_outside']['ann'], results['T_outside']['n_windows'])
    draw_quartet_tree(axs[0,1], results['T_inside']['split'], 'B. Inside 2La', results['T_inside']['ann'], results['T_inside']['n_windows'], extra_note='q(T_outside)=0.409')
    draw_quartet_tree(axs[1,0], results['T_all']['split'], 'C. All 2L windows', results['T_all']['ann'], results['T_all']['n_windows'])
    ax=axs[1,1]
    xs=[int(r['m_inside']) for r in summary_rows]
    ys=[float(r['mean_q_outside_topology']) for r in summary_rows]
    yerr=[float(r['sd_q_outside_topology']) for r in summary_rows]
    ymin=[float(r['min_q_outside_topology']) for r in summary_rows]
    ymax=[float(r['max_q_outside_topology']) for r in summary_rows]
    ax.plot(xs, ys, marker='o', color='#1f78b4', lw=2)
    ax.errorbar(xs, ys, yerr=yerr, fmt='none', ecolor='#1f78b4', alpha=0.45, capsize=2)
    ax.fill_between(xs, ymin, ymax, color='#1f78b4', alpha=0.12, linewidth=0)
    ax.set_xscale('symlog', linthresh=2)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(x) for x in xs], rotation=45, ha='right')
    ax.set_ylim(0.62,0.96)
    ax.set_xlabel('Number of 2La windows added')
    ax.set_ylabel('Quartet support for outside topology')
    ax.set_title('D. Effect of linked 2La windows', loc='left', fontsize=12, fontweight='bold')
    ax.grid(True, axis='y', color='0.9', lw=0.8)
    ax.text(0.03,0.06,'Topology remained T_outside\nfor all mixed treatments', transform=ax.transAxes, fontsize=9, bbox=dict(facecolor='white', edgecolor='0.82', pad=3))
    ax.axhline(0.408833, color='#d95f02', ls=':', lw=1.2)
    ax.text(0.98,0.13,'inside-only q(T_outside)=0.409', transform=ax.transAxes, ha='right', va='bottom', fontsize=8.5, color='#b85c00')
    ax.text(0.98,0.94,'CU: 2.245 → 0.801', transform=ax.transAxes, ha='right', va='top', fontsize=9)
    fig.suptitle('Stage 3R: species-tree sensitivity to linked 2La local genealogies', fontsize=13)
    fig.tight_layout(rect=[0,0,1,0.965])
    fig.savefig(FIG)
    fig.savefig(FIGPNG,dpi=300)
    plt.close(fig)

def write_stage3_texts(results,fixed,down_rows,classification):
    outside=topology_label(results['T_outside']['split']); inside=topology_label(results['T_inside']['split']); alltop=topology_label(results['T_all']['split'])
    fixed_by={r['treatment']:r for r in fixed}
    q_out=float(fixed_by['T_outside']['q_baseline']); q_all=float(fixed_by['T_all']['q_baseline']); q_in_base=float(fixed_by['T_inside']['q_baseline'])
    q_in_win=float(results['T_inside']['ann']['q1']); q_in_alt2=float(results['T_inside']['ann']['q2']); q_in_alt3=float(results['T_inside']['ann']['q3']); pp_in=float(results['T_inside']['ann']['pp1'])
    cu_out=float(results['T_outside']['ann']['branch_length']); cu_all=float(results['T_all']['ann']['branch_length'])
    q_reduction=(q_out-q_all)/q_out
    alt_out=1-q_out; alt_all=1-q_all
    cu_reduction=(cu_out-cu_all)/cu_out
    match_by_m={}
    for m in M_VALUES:
        vals=[as_bool(r['matches_T_outside']) for r in down_rows if str(r['m_inside'])==str(m)]
        match_by_m[str(m)]=sum(vals)/len(vals) if vals else math.nan
    METHODS.write_text(f"""# Stage 3R methods text

We tested whether the strongly arrangement-associated local genealogies within 2La influenced species-level summary-tree inference when windows were treated as local tree inputs. We reused the Stage-2R 53-tip NJ trees for usable, non-boundary 50-kb windows and mapped individuals to four biological species: arabiensis, coluzzii, gambiae, and quadriannulatus. Arrangement labels were retained for diagnostics but were not used as species in the primary ASTRAL4 run.

ASTRAL4 v1.25.4.8 was run with the documented `-a/--mapping` gene-to-species map. Treatments were all usable non-boundary 2L windows, outside-only collinear windows, and inside-only 2La windows. The outside-only topology was frozen as the empirical collinear baseline before downweighting. Downweighting retained all outside windows and added deterministic random subsets of inside windows with seeds `20261003 + 1000*m + replicate_id`; {REPS} replicates were used for finite m because full ASTRAL external-process runs were not cheap enough for 100 replicates in this environment.
""")
    RESTEXT.write_text(f"""# Stage 3R results text

The empirical collinear baseline `T_outside` was `{outside}`. The inside-only topology was `{inside}`, and the all-window topology was `{alltop}`. `T_all` matched `T_outside`, whereas `T_inside` differed from `T_outside`.

Quartet support for the outside-baseline topology declined from `{q_out:.6f}` in the outside-only analysis to `{q_all:.6f}` when all linked 2La windows were added, a relative reduction of `{100*q_reduction:.1f}%`. Combined alternative quartet support increased from `{alt_out:.4f}` outside to `{alt_all:.4f}` in the all-window analysis. The requested CU internal branch length also dropped from `{cu_out:.5f}` to `{cu_all:.6f}`, a reduction of `{100*cu_reduction:.1f}%`. CU lengths are interpreted only as summary-coalescent sensitivity metrics, not calibrated times.

Using explicit split-aware notation, `q(T_outside | outside)= {q_out:.6f}`, `q(T_outside | inside)= {q_in_base:.6f}`, and `q(T_outside | all)= {q_all:.6f}`. The inside-only analysis inferred the alternative topology `{inside}` with `q(T_inside | inside)= {q_in_win:.6f}` and near-unit ASTRAL local posterior support (`localPP={pp_in:.5f}`); the remaining third quartet had support `{q_in_alt3:.6f}`. Thus the inside-winning quartet frequency is only modestly higher than the outside-baseline quartet frequency on the same inside-only data. This illustrates sensitivity to linked-window replication: many physically linked 50-kb windows inside 2La can yield near-unit summary support even when the leading raw quartet-frequency advantage is modest.

The downweighting experiment measured whether linked 2La windows could shift species-level summary-tree inference when all outside windows were retained. The inferred topology matched `T_outside` for every mixed treatment, including all inside windows: {', '.join(f'{m}: {match_by_m[str(m)]:.2f}' for m in M_VALUES)}. The support trend is summarized in `stage3r_downweighting_summary.tsv` and plotted in the Stage 3R figure.

Stage 3R classification: **{classification}**.

This does not replace the Stage-2R primary result. Stage 2R remains the main Anopheles result: within 2La, same-arrangement samples across species are closer than same-species samples carrying opposite arrangements.
""")
    CAPTION.write_text(f"""# Stage 3R figure caption

Stage 3R species-level ASTRAL sensitivity analysis. (A) Outside-2La local trees recover the collinear-background topology `{outside}` with high quartet support and a long requested CU internal branch. (B) Inside-2La windows recover an alternative topology `{inside}`. The inside-only analysis yields near-unit ASTRAL local posterior support even though its leading quartet frequency is only modestly higher than that of the outside-baseline topology, illustrating sensitivity to linked-window replication when many physically linked 50-kb windows are treated as separate local-tree observations. (C) Combining all usable windows restores the outside topology but reduces its quartet support and requested CU branch length. (D) Quartet support is always defined with respect to the fixed outside-2La topology, irrespective of which split is ranked first in an individual ASTRAL run; this support declines as increasing numbers of linked 2La windows are added while the inferred mixed-treatment topology remains unchanged. The 429-window endpoint represents all outside windows plus all usable 2La windows, whereas the inside-only analysis is shown separately in panel B. CU length is a summary-coalescent sensitivity metric here and is not calibrated divergence time.
""")
    REPORT.write_text(f"""# Stage 3R — species-level summary-tree sensitivity to 2La-linked local genealogies

Stage 3R asked whether hundreds of linked arrangement-associated local genealogies inside 2La propagate into species-level summary-tree inference relative to the collinear 2L background.

- `T_outside`: `{outside}`
- `T_inside`: `{inside}`
- `T_all`: `{alltop}`
- final classification: **{classification}**

`T_inside` differs from `T_outside`, but `T_all` equals `T_outside`. Adding all 2La windows reduces `q(T_outside)` from approximately `{q_out:.3f}` to `{q_all:.3f}` and reduces the requested CU internal branch length from approximately `{cu_out:.3f}` to `{cu_all:.3f}`. The full combined topology is therefore stable, but the linked inversion windows weaken the support/branch-length profile. Inside-only linked windows nevertheless yield near-unit support for an alternative topology (`localPP={pp_in:.5f}`), with `q(T_inside | inside)={q_in_win:.3f}` versus `q(T_outside | inside)={q_in_base:.3f}`.

Interpretation: the extreme local arrangement signal does not overturn the four-species all-window summary topology, but it materially affects summary support and CU branch-length sensitivity. The Stage-2R crossed MalariaGEN result remains the primary Anopheles biological result.

ASTRAL branch lengths, when present, are treated only as summary-coalescent sensitivity metrics. SU lengths are not interpreted as calibrated substitution lengths or divergence times.
""")

def update_docs(classification,results,down_rows):
    block=f"""

## Stage 3R — ASTRAL species-tree sensitivity

Stage 3R reused the Stage-2R local 53-individual NJ trees and ASTRAL4's documented multi-individual mapping (`-a`) to test whether linked arrangement-dominated 2La windows alter four-species summary-tree inference. The Stage-2R phased validation text was corrected to 6/6 concordant deterministic validation windows, and the main Anopheles figure was rebuilt with matplotlib axes.

- `T_outside`: {topology_label(results['T_outside']['split'])}
- `T_inside`: {topology_label(results['T_inside']['split'])}
- `T_all`: {topology_label(results['T_all']['split'])}
- classification: **{classification}**

The primary Anopheles biological result remains Stage 2R's crossed MalariaGEN result inside 2La; Stage 3R is a sensitivity analysis of propagation into species-level summary-tree inference. The improved Stage 3R figure shows that 2La does not flip the combined species topology, but progressively erodes its quartet and CU support, while 2La-only windows support a different topology with near-unit ASTRAL support.
"""
    for p in [README,PROJECT_STATUS]:
        old=p.read_text() if p.exists() else ''
        marker='## Stage 3R — ASTRAL species-tree sensitivity'
        if marker in old: old=old.split(marker)[0].rstrip()+'\n'
        p.write_text(old.rstrip()+block)


def update_existing_outputs():
    fix_stage2_texts(); check_stage2_numbers(); make_main_figure()
    results=load_existing_astral_results()
    fixed=read_tsv(FIXED)
    down_rows=read_tsv(DOWN)
    classification='SUPPORT/BRANCH-LENGTH EFFECT'
    make_stage3_fig(results, down_rows, fixed)
    write_stage3_texts(results, fixed, down_rows, classification)
    update_docs(classification, results, down_rows)
    outputs=[DOWN_SUMMARY,RESTEXT,CAPTION,REPORT,FIG,FIGPNG,MAIN_PDF,MAIN_PNG,README,PROJECT_STATUS]
    old=json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
    old.update({'stage':'Stage 3R','classification':classification,'T_outside':topology_label(results['T_outside']['split']),'T_inside':topology_label(results['T_inside']['split']),'T_all':topology_label(results['T_all']['split']),'q_Toutside_outside':0.931104,'q_Toutside_inside':0.408833,'q_Tinside_inside':0.456944,'q_Toutside_all':0.701540,'support_definition':'canonical fixed T_outside split: arabiensis,quadriannulatus|coluzzii,gambiae','outputs':{**old.get('outputs',{}), **{str(p.relative_to(REPO)):sha256(p) for p in outputs if p.exists()}}})
    MANIFEST.write_text(json.dumps(old,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'updated_existing':True,'q_Toutside_outside':0.931104,'q_Toutside_inside':0.408833,'q_Tinside_inside':0.456944,'q_Toutside_all':0.701540},sort_keys=True))
    return 0

def run(args):
    fix_stage2_texts(); check_stage2_numbers(); make_main_figure()
    grid,qc,inside,outside,clean=load_stage2(); trees=load_trees(); mapping=write_species_map(); audit_trees(trees,clean,inside,outside,mapping)
    results,fixed,mapfile,infiles=do_astral(trees,inside,outside,clean,mapping)
    down_rows=downweight(trees,inside,outside,results,mapfile)
    arrangement_control()
    changed_all=results['T_all']['split']!=results['T_outside']['split']; changed_inside=results['T_inside']['split']!=results['T_outside']['split']
    classification='SPECIES-TOPOLOGY EFFECT' if changed_all else 'SUPPORT/BRANCH-LENGTH EFFECT' if changed_inside else 'LOCAL-ONLY EFFECT'
    make_stage3_fig(results,down_rows); write_stage3_texts(results,fixed,down_rows,classification); update_docs(classification,results,down_rows)
    outputs=[MAP,AUDIT,ENV,TOPO_SUM,FIXED,DOWN,DOWN_SUMMARY,BRANCH,ARRCTRL,METHODS,RESTEXT,CAPTION,REPORT,FIG,FIGPNG,MAIN_PDF,MAIN_PNG,STAGE2_CAPTION,STAGE2_RESULTS,README,PROJECT_STATUS]
    for p in ASTRAL_DIR.glob('*'): outputs.append(p)
    MANIFEST.write_text(json.dumps({'stage':'Stage 3R','classification':classification,'astral4':str(ASTRAL4),'direct_mapping_supported':True,'n_inside_trees':len(inside),'n_outside_trees':len(outside),'n_clean_trees':len(clean),'T_outside':topology_label(results['T_outside']['split']),'T_inside':topology_label(results['T_inside']['split']),'T_all':topology_label(results['T_all']['split']),'stage2_haplotype_concordance':'6/6','stage2_core_numbers_unchanged':True,'downweight_reps_finite_m':REPS,'outputs':{str(p.relative_to(REPO)):sha256(p) for p in outputs if p.exists()}},indent=2,sort_keys=True)+'\n')
    print(json.dumps({'classification':classification,'T_outside':topology_label(results['T_outside']['split']),'T_inside':topology_label(results['T_inside']['split']),'T_all':topology_label(results['T_all']['split']),'inside':len(inside),'outside':len(outside)},sort_keys=True))
    return 0

class Tests(unittest.TestCase):
    def test_haplotype_6_of_6(self):
        rows=read_tsv(HAPVAL); self.assertEqual(len(rows),6); self.assertEqual(sum(as_bool(r['concordant_direction']) for r in rows),6)
    def test_stage2_numbers(self): check_stage2_numbers()
    def test_species_map_counts(self):
        rows=[r for r in read_tsv(PRIMARY) if r['include_primary']=='true']; self.assertEqual(Counter(r['species'] for r in rows),Counter({'arabiensis':12,'coluzzii':11,'gambiae':20,'quadriannulatus':10}))
    def test_clean_counts(self):
        _,_,inside,outside,clean=load_stage2(); self.assertEqual(len(inside),429); self.assertEqual(len(outside),547); self.assertEqual(len(clean),976)

    def test_downweighting_summary_endpoints(self):
        if DOWN_SUMMARY.exists():
            rows={r['m_inside']:r for r in read_tsv(DOWN_SUMMARY)}
            self.assertAlmostEqual(float(rows['0']['mean_q_outside_topology']),0.931104,places=6)
            self.assertAlmostEqual(float(rows['429']['mean_q_outside_topology']),0.70154,places=5)
            self.assertNotAlmostEqual(float(rows['429']['mean_q_outside_topology']),0.408833,places=3)
            self.assertEqual(rows['0']['support_source'],'T_outside exact')
            self.assertEqual(rows['429']['support_source'],'T_all exact')
            self.assertTrue(all(abs(float(r['fraction_topology_matches_outside'])-1.0)<1e-12 for r in rows.values()))
    def test_canonical_baseline_support(self):
        results=load_existing_astral_results()
        fixed=fixed_score_rows_by_treatment()
        self.assertEqual(results['T_inside']['split'], T_INSIDE_SPLIT)
        self.assertNotEqual(results['T_inside']['split'], T_OUTSIDE_SPLIT)
        self.assertAlmostEqual(float(results['T_inside']['ann']['q1']),0.456944,places=6)
        self.assertAlmostEqual(baseline_quartet_support('T_inside', results['T_inside']['split'], results['T_inside']['ann'], fixed)[0],0.408833,places=6)

def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument('--run-tests',action='store_true'); ap.add_argument('--update-existing', action='store_true', help='Regenerate Stage 3R figure/prose from existing ASTRAL outputs without rerunning ASTRAL.'); args=ap.parse_args(argv)
    if args.run_tests:
        res=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests)); return 0 if res.wasSuccessful() else 1
    if args.update_existing:
        return update_existing_outputs()
    return run(args)
if __name__=='__main__': raise SystemExit(main())
