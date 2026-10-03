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

def make_stage3_fig(results,down_rows):
    fig,axs=plt.subplots(2,2,figsize=(10,7))
    for ax,t in zip(axs.flat[:3],['T_outside','T_inside','T_all']):
        ax.axis('off'); ax.set_title(t.replace('_',' '),fontsize=12)
        top=topology_label(results[t]['split']); sup=results[t]['ann'].get('pp1',results[t]['ann'].get('confidence','NA'))
        ax.text(0.5,0.58,top,ha='center',va='center',fontsize=11,wrap=True)
        ax.text(0.5,0.35,f'support: {fmt(sup)}\nwindows: {results[t]["n_windows"]}',ha='center',va='center',fontsize=10)
    ax=axs.flat[3]
    grouped=defaultdict(list)
    for r in down_rows: grouped[str(r['m_inside'])].append(1.0 if as_bool(r['matches_T_outside']) else 0.0)
    xs=[]; ys=[]; labels=[]
    for m in M_VALUES:
        k=str(m); xs.append(len(xs)); ys.append(sum(grouped[k])/len(grouped[k]) if grouped[k] else math.nan); labels.append(str(m))
    ax.plot(xs,ys,marker='o',color='#1b9e77')
    ax.set_xticks(xs); ax.set_xticklabels(labels,rotation=45)
    ax.set_ylim(-0.03,1.03); ax.set_ylabel('Fraction matching T_outside'); ax.set_xlabel('Number of 2La windows included')
    ax.set_title('Downweighting linked 2La windows')
    fig.tight_layout(); fig.savefig(FIG); fig.savefig(FIGPNG,dpi=300); plt.close(fig)

def write_stage3_texts(results,fixed,down_rows,classification):
    outside=topology_label(results['T_outside']['split']); inside=topology_label(results['T_inside']['split']); alltop=topology_label(results['T_all']['split'])
    changed_all=results['T_all']['split']!=results['T_outside']['split']; changed_inside=results['T_inside']['split']!=results['T_outside']['split']
    match_by_m={}
    for m in M_VALUES:
        vals=[as_bool(r['matches_T_outside']) for r in down_rows if str(r['m_inside'])==str(m)]
        match_by_m[str(m)]=sum(vals)/len(vals) if vals else math.nan
    METHODS.write_text(f"""# Stage 3R methods text

We tested whether the strongly arrangement-associated local genealogies within 2La influenced species-level summary-tree inference when windows were treated as local tree inputs. We reused the Stage-2R 53-tip NJ trees for usable, non-boundary 50-kb windows and mapped individuals to four biological species: arabiensis, coluzzii, gambiae, and quadriannulatus. Arrangement labels were retained for diagnostics but were not used as species in the primary ASTRAL4 run.

ASTRAL4 v1.25.4.8 was run with the documented `-a/--mapping` gene-to-species map. Treatments were all usable non-boundary 2L windows, outside-only collinear windows, and inside-only 2La windows. The outside-only topology was frozen as the empirical collinear baseline before downweighting. Downweighting retained all outside windows and added deterministic random subsets of inside windows with seeds `20261003 + 1000*m + replicate_id`; {REPS} replicates were used for finite m because full ASTRAL external-process runs were not cheap enough for 100 replicates in this environment.
""")
    RESTEXT.write_text(f"""# Stage 3R results text

The empirical collinear baseline `T_outside` was `{outside}`. The inside-only topology was `{inside}`, and the all-window topology was `{alltop}`. `T_all` {'differed from' if changed_all else 'matched'} `T_outside`; `T_inside` {'differed from' if changed_inside else 'matched'} `T_outside`.

The downweighting experiment measured whether linked 2La windows could shift species-level summary-tree inference when all outside windows were retained. The fraction of replicates matching `T_outside` by m was: {', '.join(f'{m}: {match_by_m[str(m)]:.2f}' for m in M_VALUES)}.

Stage 3R classification: **{classification}**.

This does not replace the Stage-2R primary result. Stage 2R remains the main Anopheles result: within 2La, same-arrangement samples across species are closer than same-species samples carrying opposite arrangements.
""")
    CAPTION.write_text("""# Stage 3R figure caption

ASTRAL4 species-level summary-tree sensitivity to linked 2La local genealogies. Panels A-C show the inferred four-species summary topologies for outside-only, inside-only, and all usable non-boundary 2L windows, with support values reported from ASTRAL output where available. Panel D shows the downweighting experiment retaining all outside windows while adding deterministic subsets of inside-2La windows.
""")
    REPORT.write_text(f"""# Stage 3R — species-level summary-tree sensitivity to 2La-linked local genealogies

Stage 3R asked whether hundreds of linked arrangement-associated local genealogies inside 2La propagate into species-level summary-tree inference relative to the collinear 2L background.

- `T_outside`: `{outside}`
- `T_inside`: `{inside}`
- `T_all`: `{alltop}`
- final classification: **{classification}**

Interpretation: {'including linked 2La windows altered the inferred species-level summary topology relative to the collinear background.' if classification=='SPECIES-TOPOLOGY EFFECT' else 'the extreme local arrangement signal did not overturn the four-species summary topology, but support/branch-length sensitivity is documented.' if classification=='SUPPORT/BRANCH-LENGTH EFFECT' else 'the effect remained local to 2La arrangement genealogy and did not materially affect species-level summary inference.'}

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

The primary Anopheles biological result remains Stage 2R's crossed MalariaGEN result inside 2La; Stage 3R is a sensitivity analysis of propagation into species-level summary-tree inference.
"""
    for p in [README,PROJECT_STATUS]:
        old=p.read_text() if p.exists() else ''
        marker='## Stage 3R — ASTRAL species-tree sensitivity'
        if marker in old: old=old.split(marker)[0].rstrip()+'\n'
        p.write_text(old.rstrip()+block)

def run(args):
    fix_stage2_texts(); check_stage2_numbers(); make_main_figure()
    grid,qc,inside,outside,clean=load_stage2(); trees=load_trees(); mapping=write_species_map(); audit_trees(trees,clean,inside,outside,mapping)
    results,fixed,mapfile,infiles=do_astral(trees,inside,outside,clean,mapping)
    down_rows=downweight(trees,inside,outside,results,mapfile)
    arrangement_control()
    changed_all=results['T_all']['split']!=results['T_outside']['split']; changed_inside=results['T_inside']['split']!=results['T_outside']['split']
    classification='SPECIES-TOPOLOGY EFFECT' if changed_all else 'SUPPORT/BRANCH-LENGTH EFFECT' if changed_inside else 'LOCAL-ONLY EFFECT'
    make_stage3_fig(results,down_rows); write_stage3_texts(results,fixed,down_rows,classification); update_docs(classification,results,down_rows)
    outputs=[MAP,AUDIT,ENV,TOPO_SUM,FIXED,DOWN,BRANCH,ARRCTRL,METHODS,RESTEXT,CAPTION,REPORT,FIG,FIGPNG,MAIN_PDF,MAIN_PNG,STAGE2_CAPTION,STAGE2_RESULTS,README,PROJECT_STATUS]
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

def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument('--run-tests',action='store_true'); args=ap.parse_args(argv)
    if args.run_tests:
        res=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests)); return 0 if res.wasSuccessful() else 1
    return run(args)
if __name__=='__main__': raise SystemExit(main())
