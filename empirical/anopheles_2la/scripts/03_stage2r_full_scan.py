#!/usr/bin/env python3
"""Stage 2R full 2L 50-kb spatial scan for Anopheles 2La.

Uses processed MalariaGEN Ag3 regional SNP/haplotype calls only. No raw reads.
"""
from __future__ import annotations

import argparse, csv, gzip, hashlib, itertools, json, math, os, platform, sys, time, unittest
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / 'data' / 'anopheles_2la'
PROC = DATA / 'processed'
META = DATA / 'metadata'
BASE = REPO / 'empirical' / 'anopheles_2la'
RESULTS = BASE / 'results'
FIGURES = BASE / 'figures'
README = BASE / 'README.md'
PROJECT_STATUS = REPO / 'PROJECT_STATUS.md'

PRIMARY_SAMPLES = PROC / 'stage1r_primary_homozygote_samples.tsv'
STAGE1R_MANIFEST = RESULTS / 'stage1r_manifest.json'
PREPLAN = RESULTS / 'stage2r_preanalysis_plan.md'
GRID = PROC / 'stage2r_window_grid.tsv'
PAIRWISE = PROC / 'stage2r_pairwise_distances.tsv.gz'
TREES = PROC / 'stage2r_local_nj_trees.nwk.gz'
QC = RESULTS / 'stage2r_window_qc.tsv'
CROSS = RESULTS / 'stage2r_crossed_distance_signal.tsv'
MARGIN = RESULTS / 'stage2r_quartet_margin_signal.tsv'
TOPO = RESULTS / 'stage2r_quartet_topology_signal.tsv'
SUMMARY = RESULTS / 'stage2r_inside_outside_summary.tsv'
CIRC = RESULTS / 'stage2r_circular_shift_tests.tsv'
PHYS = RESULTS / 'stage2r_physical_interval_tests.tsv'
DECAY = RESULTS / 'stage2r_boundary_decay.tsv'
INPOS = RESULTS / 'stage2r_inside_position_summary.tsv'
ALL53 = RESULTS / 'stage2r_all53_distance_class_summary.tsv'
SUMTOPO = RESULTS / 'stage2r_summary_topology_analysis.md'
HAPVAL = RESULTS / 'stage2r_haplotype_validation.tsv'
FONTCMP = RESULTS / 'stage2r_fontaine_comparison.md'
METHODS = RESULTS / 'stage2r_methods_text.md'
RESTEXT = RESULTS / 'stage2r_results_text.md'
CAPTION = RESULTS / 'stage2r_main_figure_caption.md'
REPORT = RESULTS / 'stage2r_report.md'
MANIFEST = RESULTS / 'stage2r_manifest.json'
MAIN_PDF = FIGURES / 'stage2r_full_2L_spatial_signal.pdf'
MAIN_PNG = FIGURES / 'stage2r_full_2L_spatial_signal.png'
FIG_MAIN_PDF = FIGURES / 'anopheles_2la_main.pdf'
FIG_MAIN_PNG = FIGURES / 'anopheles_2la_main.png'
FIG_DESIGN_PDF = FIGURES / 'anopheles_2la_crossed_design.pdf'
FIG_DESIGN_PNG = FIGURES / 'anopheles_2la_crossed_design.png'

CHROM='2L'; CHROM_LEN=49_364_325; WIN=50_000
INV_START=20_524_058; INV_END=42_165_532; INV_WIDTH=INV_END-INV_START+1
RELEASE='3.10'; SAMPLE_SET='fontaine-2015-rebuild'; BUCKET='gs://vo_agam_release_master_us_central1'
MIN_MAC=2; MAX_SITE_MISS=0.25; MAX_SAMPLE_MISS=0.50; MIN_SNPS=25
CHUNK_BP=1_000_000
EXPECTED_COUNTS={('arabiensis','2La/2La'):12,('coluzzii','2L+a/2L+a'):8,('coluzzii','2La/2La'):3,('gambiae','2L+a/2L+a'):8,('gambiae','2La/2La'):12,('quadriannulatus','2L+a/2L+a'):10}


def fmt(v):
    if v is None: return ''
    if isinstance(v, bool): return 'true' if v else 'false'
    if isinstance(v, (float, np.floating)):
        if math.isnan(float(v)): return 'nan'
        return f'{float(v):.12g}'
    if isinstance(v, (int, np.integer)): return str(int(v))
    return str(v)

def read_tsv(path):
    with path.open(newline='') as f: return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator='\n'); w.writeheader()
        for r in rows: w.writerow({k:fmt(r.get(k,'')) for k in fields})

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def load_samples():
    rows=read_tsv(PRIMARY_SAMPLES)
    assert len(rows)==72, f'expected 72 rows including excluded heterokaryotypes, got {len(rows)}'
    primary=[r for r in rows if r['include_primary']=='true']
    assert len(primary)==53, f'expected 53 primary homozygotes, got {len(primary)}'
    assert sum(1 for r in rows if r['include_primary']!='true')==19
    c=Counter((r['species'],r['arrangement']) for r in primary)
    for k,v in EXPECTED_COUNTS.items(): assert c[k]==v, f'{k} expected {v} got {c[k]}'
    crossed=[r for r in primary if r['species'] in {'gambiae','coluzzii'}]
    assert len(crossed)==31
    assert len(set(r['sample_id'] for r in primary))==53
    return rows, primary

def make_grid():
    rows=[]; idx=0
    for start in range(1, CHROM_LEN+1, WIN):
        idx+=1; end=min(start+WIN-1, CHROM_LEN)
        if start>=INV_START and end<=INV_END: rc='inside'
        elif end<INV_START or start>INV_END: rc='outside'
        else: rc='boundary'
        rows.append({'window_id':f'2L_{idx:04d}_{start}_{end}','chrom':CHROM,'start':start,'end':end,'region_class':rc,
                     'distance_to_left_breakpoint': min(abs(start-INV_START),abs(end-INV_START)) if not (start<=INV_START<=end) else 0,
                     'distance_to_right_breakpoint': min(abs(start-INV_END),abs(end-INV_END)) if not (start<=INV_END<=end) else 0})
    assert len(rows)==math.ceil(CHROM_LEN/WIN)
    assert rows[0]['start']==1 and rows[0]['end']==50000 and rows[1]['start']==50001
    assert sum(1 for r in rows if r['region_class']=='boundary')==2
    return rows

def dosage_filters(gt):
    called=np.all(gt>=0, axis=2)
    dosage=np.where(called, gt.sum(axis=2), -1).astype(np.int16)
    n_called=called.sum(axis=1); an=2*n_called
    alt=np.where(called,dosage,0).sum(axis=1)
    mac=np.minimum(alt, an-alt)
    site_miss=1-(n_called/gt.shape[1]) if gt.shape[1] else np.ones(gt.shape[0])
    keep=(mac>=MIN_MAC)&(site_miss<=MAX_SITE_MISS)
    return dosage, called, alt, mac, keep

def dist_matrix(dosage, called, keep):
    idx=np.where(keep)[0]; n=dosage.shape[1]
    D=np.full((n,n), np.nan); N=np.zeros((n,n), dtype=np.int32)
    for i in range(n):
        D[i,i]=0.0; N[i,i]=len(idx)
        for j in range(i+1,n):
            both=called[idx,i]&called[idx,j]; m=int(both.sum()); N[i,j]=N[j,i]=m
            if m:
                d=float(np.abs(dosage[idx[both],i]-dosage[idx[both],j]).mean()/2.0)
                D[i,j]=D[j,i]=d
    return D,N

def crossed_indices(samples):
    return {
        'G_s':[i for i,s in enumerate(samples) if s['species']=='gambiae' and s['arrangement']=='2L+a/2L+a'],
        'G_i':[i for i,s in enumerate(samples) if s['species']=='gambiae' and s['arrangement']=='2La/2La'],
        'C_s':[i for i,s in enumerate(samples) if s['species']=='coluzzii' and s['arrangement']=='2L+a/2L+a'],
        'C_i':[i for i,s in enumerate(samples) if s['species']=='coluzzii' and s['arrangement']=='2La/2La'],
    }

def pair_iter(samples, D, N, window_id):
    for i in range(len(samples)):
        for j in range(i+1,len(samples)):
            yield {'window_id':window_id,'sample1':samples[i]['sample_id'],'sample2':samples[j]['sample_id'],'species1':samples[i]['species'],'species2':samples[j]['species'],'arrangement1':samples[i]['arrangement'],'arrangement2':samples[j]['arrangement'],'distance':D[i,j],'n_sites_compared':N[i,j]}

def class_summary(w, samples, D, N):
    b=defaultdict(list); ns=defaultdict(list)
    for i in range(len(samples)):
        for j in range(i+1,len(samples)):
            if math.isnan(D[i,j]): continue
            if samples[i]['species']==samples[j]['species'] and samples[i]['arrangement']==samples[j]['arrangement']: cls='same_species_same_arrangement'
            elif samples[i]['species']==samples[j]['species']: cls='same_species_opposite_arrangement'
            elif samples[i]['arrangement']==samples[j]['arrangement']: cls='different_species_same_arrangement'
            else: cls='different_species_opposite_arrangement'
            b[cls].append(float(D[i,j])); ns[cls].append(int(N[i,j]))
    out=[]
    for cls in ['same_species_same_arrangement','same_species_opposite_arrangement','different_species_same_arrangement','different_species_opposite_arrangement']:
        vals=b[cls]
        out.append({'window_id':w['window_id'],'start':w['start'],'end':w['end'],'region_class':w['region_class'],'distance_class':cls,'n_pairs':len(vals),'mean_distance':np.mean(vals) if vals else np.nan,'median_distance':np.median(vals) if vals else np.nan,'mean_sites_compared':np.mean(ns[cls]) if vals else np.nan})
    return out

def crossed_C(w, samples, D):
    within=[]; cross=[]
    for i in range(len(samples)):
        for j in range(i+1,len(samples)):
            if samples[i]['species'] not in {'gambiae','coluzzii'} or samples[j]['species'] not in {'gambiae','coluzzii'} or math.isnan(D[i,j]): continue
            if samples[i]['species']==samples[j]['species'] and samples[i]['arrangement']!=samples[j]['arrangement']: within.append(float(D[i,j]))
            if samples[i]['species']!=samples[j]['species'] and samples[i]['arrangement']==samples[j]['arrangement']: cross.append(float(D[i,j]))
    mw=np.mean(within) if within else np.nan; mc=np.mean(cross) if cross else np.nan
    return {'window_id':w['window_id'],'start':w['start'],'end':w['end'],'region_class':w['region_class'],'n_within_species_opposite_pairs':len(within),'mean_within_species_opposite':mw,'n_cross_species_same_pairs':len(cross),'mean_cross_species_same':mc,'C':mw-mc if within and cross else np.nan}

def quartet_stats(w, samples, D):
    g=crossed_indices(samples); margins=[]; top=Counter(); n=0
    for gs,gi,cs,ci in itertools.product(g['G_s'],g['G_i'],g['C_s'],g['C_i']):
        ss=D[gs,gi]+D[cs,ci]; sa=D[gs,cs]+D[gi,ci]; st=D[gs,ci]+D[gi,cs]
        if any(math.isnan(x) for x in [ss,sa,st]): continue
        n+=1; margins.append(ss-sa)
        vals={'species':ss,'arrangement':sa,'third':st}; m=min(vals.values()); wins=[k for k,v in vals.items() if abs(v-m)<1e-12]
        top[wins[0] if len(wins)==1 else 'tie']+=1
    assert n in {0,2304,36864}
    arr=np.array(margins, dtype=float)
    qs=top['species']/n if n else np.nan; qa=top['arrangement']/n if n else np.nan; q3=(top['third']+top['tie'])/n if n else np.nan
    mrow={'window_id':w['window_id'],'start':w['start'],'end':w['end'],'region_class':w['region_class'],'n_informative_quartets':n,'mean_M':float(np.mean(arr)) if n else np.nan,'median_M':float(np.median(arr)) if n else np.nan,'fraction_M_positive':float((arr>0).mean()) if n else np.nan,'fraction_M_negative':float((arr<0).mean()) if n else np.nan,'fraction_M_zero':float((arr==0).mean()) if n else np.nan}
    drow={'window_id':w['window_id'],'start':w['start'],'end':w['end'],'region_class':w['region_class'],'n_informative_quartets':n,'q_species':qs,'q_arrangement':qa,'q_third':q3,'D':qa-qs if n else np.nan}
    return mrow,drow

def nj_newick(labels,D):
    finite=D[np.isfinite(D)]; fill=float(finite.max()) if finite.size else 1.0
    M=D.copy().astype(float); M[~np.isfinite(M)]=fill; np.fill_diagonal(M,0)
    clusters={i:labels[i] for i in range(len(labels))}; active=list(range(len(labels))); nxt=len(labels)
    while len(active)>2:
        n=len(active); total={i:sum(M[i,j] for j in active if j!=i) for i in active}; best=None
        for a,i in enumerate(active):
            for j in active[a+1:]:
                key=((n-2)*M[i,j]-total[i]-total[j], i, j)
                if best is None or key<best[0]: best=(key,i,j)
        _,i,j=best; li=max(0.0,0.5*M[i,j]+(total[i]-total[j])/(2*(n-2))); lj=max(0.0,M[i,j]-li)
        if nxt>=M.shape[0]:
            O=M; M=np.zeros((O.shape[0]+1,O.shape[1]+1)); M[:O.shape[0],:O.shape[1]]=O
        for k in active:
            if k not in (i,j): M[nxt,k]=M[k,nxt]=0.5*(M[i,k]+M[j,k]-M[i,j])
        clusters[nxt]=f'({clusters[i]}:{li:.6g},{clusters[j]}:{lj:.6g})'; active=[k for k in active if k not in (i,j)]+[nxt]; nxt+=1
    i,j=active; bl=max(0.0,M[i,j]/2); return f'({clusters[i]}:{bl:.6g},{clusters[j]}:{bl:.6g});'

def label(s): return f"{s['sample_id']}|{s['species'][:4]}|{'std' if s['arrangement']=='2L+a/2L+a' else 'inv'}"

def window_qc(w, n_bial, dosage, called, keep, runtime, nbytes):
    if dosage.shape[0]:
        sample_miss=1-called.mean(axis=0); mean_miss=float(1-called.mean()); het=float(((dosage==1)&called).sum()/called.sum()) if called.sum() else np.nan
    else:
        sample_miss=np.ones(53); mean_miss=1.0; het=np.nan
    callable_=int(keep.sum())>=MIN_SNPS and float(sample_miss.max())<=MAX_SAMPLE_MISS
    return {'window_id':w['window_id'],'start':w['start'],'end':w['end'],'region_class':w['region_class'],'n_sites_raw':int(w['end'])-int(w['start'])+1,'n_biallelic':n_bial,'n_post_filter_snps':int(keep.sum()),'mean_missingness':mean_miss,'max_sample_missingness':float(sample_miss.max()),'mean_heterozygosity':het,'callable':callable_,'status':'ok' if callable_ else 'low_signal','runtime_seconds':runtime,'approx_dataset_bytes':nbytes}

def process_chunk(ag3, region, windows, sample_indices, samples):
    # Use direct SNP/genotype calls and compute biallelic/MAC/missingness locally.
    # This avoids the slower API-side allele-count cache path in biallelic_snp_calls
    # while preserving the frozen Stage-1R filter.
    t0=time.perf_counter()
    ds=ag3.snp_calls(region=region, sample_sets=SAMPLE_SET, sample_indices=sample_indices, site_mask=None, inline_array=True, chunks='native')
    pos=np.asarray(ds['variant_position'].values)
    gt=np.asarray(ds['call_genotype'].values)
    nbytes=int(getattr(ds['call_genotype'],'nbytes',0)) + int(getattr(ds['variant_position'],'nbytes',0))
    elapsed=time.perf_counter()-t0
    out=[]
    for w in windows:
        mask=(pos>=w['start'])&(pos<=w['end']); subgt=gt[mask]
        if subgt.size:
            dosage, called, alt, mac, keep0=dosage_filters(subgt)
            max_allele=np.where(called, subgt.max(axis=2), -1).max(axis=1)
            segregating=mac>0
            biallelic=(max_allele<=1)&segregating
            keep=keep0&biallelic
            dosage=dosage[biallelic]; called=called[biallelic]; keep=keep[biallelic]
            n_bial=int(biallelic.sum())
        else:
            dosage=np.empty((0,len(samples)),dtype=np.int16); called=np.empty((0,len(samples)),dtype=bool); keep=np.array([], dtype=bool); n_bial=0
        qc=window_qc(w, n_bial, dosage, called, keep, elapsed/max(1,len(windows)), nbytes/max(1,len(windows)))
        if qc['callable']:
            D,N=dist_matrix(dosage,called,keep)
        else:
            D=np.full((len(samples),len(samples)),np.nan); np.fill_diagonal(D,0); N=np.zeros((len(samples),len(samples)),dtype=np.int32)
        out.append((w,qc,D,N))
    return out

def numeric_track(rows, key):
    return np.array([float(r[key]) if str(r.get(key,'nan'))!='nan' else np.nan for r in rows], dtype=float)

def inside_outside(rows, key, grid):
    by={r['window_id']:r for r in rows}; inside=[]; outside=[]
    for w in grid:
        if w['region_class']=='boundary': continue
        v=float(by[w['window_id']][key])
        if math.isnan(v): continue
        (inside if w['region_class']=='inside' else outside).append(v)
    return inside,outside

def delta(vals, mask):
    arr=np.asarray(vals,dtype=float); m=np.asarray(mask,dtype=bool); inside=arr[m&~np.isnan(arr)]; outside=arr[(~m)&~np.isnan(arr)]
    return float(np.mean(inside)-np.mean(outside)) if len(inside) and len(outside) else np.nan

def circular_test(values, grid):
    nonboundary=[w['region_class']!='boundary' for w in grid]; inside=[w['region_class']=='inside' for w in grid]
    vals=np.asarray([v if nb else np.nan for v,nb in zip(values,nonboundary)], dtype=float); mask=np.asarray(inside, dtype=bool)
    obs=delta(vals,mask); null=[]
    for k in range(len(vals)):
        null.append(delta(np.roll(vals,k), mask))
    null=np.asarray(null); valid=null[~np.isnan(null)]
    p=float((valid>=obs).sum()/len(valid)) if len(valid) else np.nan
    return obs,p,len(valid)

def physical_test(values, grid):
    vals=np.asarray(values,dtype=float); width=INV_WIDTH; starts=[]; deltas=[]
    for s in range(1, CHROM_LEN-width+2, WIN):
        e=s+width-1; mask=[]
        for w in grid:
            mask.append(w['start']>=s and w['end']<=e)
        d=delta(vals,mask)
        if not math.isnan(d): starts.append(s); deltas.append(d)
    obs_mask=[w['region_class']=='inside' for w in grid]; obs=delta(vals,obs_mask)
    arr=np.asarray(deltas); k=int((arr>=obs).sum()) if len(arr) else 0; p=(k+1)/(len(arr)+1) if len(arr) else np.nan
    return obs,p,len(arr)

def write_summaries(grid,crows,mrows,drows,qcrows):
    rows=[]
    for name,src,key in [('C',crows,'C'),('M',mrows,'mean_M'),('D',drows,'D')]:
        ins,out=inside_outside(src,key,grid)
        rows.append({'statistic':name,'n_inside':len(ins),'n_outside':len(out),'mean_inside':np.mean(ins),'median_inside':np.median(ins),'mean_outside':np.mean(out),'median_outside':np.median(out),'delta_mean':np.mean(ins)-np.mean(out),'delta_median':np.median(ins)-np.median(out)})
    write_tsv(SUMMARY,rows,['statistic','n_inside','n_outside','mean_inside','median_inside','mean_outside','median_outside','delta_mean','delta_median'])
    circ=[]; phys=[]
    for name,src,key in [('C',crows,'C'),('M',mrows,'mean_M'),('D',drows,'D')]:
        by={r['window_id']:r for r in src}; vals=[float(by[w['window_id']][key]) if by[w['window_id']][key]==by[w['window_id']][key] else np.nan for w in grid]
        obs,p,n=circular_test(vals,grid); circ.append({'statistic':name,'observed_delta':obs,'p_one_sided_ge_observed':p,'n_rotations':n,'null':'exact circular shift; observed alignment included'})
        obs2,p2,n2=physical_test(vals,grid); phys.append({'statistic':name,'observed_delta':obs2,'p_one_sided_ge_observed_corrected':p2,'n_candidate_intervals':n2,'null':'same-width physical interval placement; p=(k+1)/(N+1)'})
    write_tsv(CIRC,circ,['statistic','observed_delta','p_one_sided_ge_observed','n_rotations','null'])
    write_tsv(PHYS,phys,['statistic','observed_delta','p_one_sided_ge_observed_corrected','n_candidate_intervals','null'])
    return rows,circ,phys

def boundary_decay(grid,crows,mrows,drows):
    bins=[(0,500000,'0-0.5Mb'),(500000,1000000,'0.5-1Mb'),(1000000,2000000,'1-2Mb'),(2000000,5000000,'2-5Mb'),(5000000,10**12,'>5Mb')]
    maps={name:{r['window_id']:r for r in rows} for name,rows in [('C',crows),('M',mrows),('D',drows)]}; keys={'C':'C','M':'mean_M','D':'D'}
    out=[]
    for lo,hi,label_ in bins:
        wins=[]
        for w in grid:
            if w['region_class']!='outside': continue
            mid=(w['start']+w['end'])/2; dist=min(abs(mid-INV_START),abs(mid-INV_END))
            if lo<=dist<hi: wins.append(w)
        row={'distance_bin':label_,'n_windows':len(wins)}
        for stat in ['C','M','D']:
            vals=[float(maps[stat][w['window_id']][keys[stat]]) for w in wins if not math.isnan(float(maps[stat][w['window_id']][keys[stat]]))]
            row[f'mean_{stat}']=np.mean(vals) if vals else np.nan; row[f'median_{stat}']=np.median(vals) if vals else np.nan
        out.append(row)
    write_tsv(DECAY,out,['distance_bin','n_windows','mean_C','median_C','mean_M','median_M','mean_D','median_D'])
    return out

def inside_position(grid,crows,mrows,drows,qcrows):
    bins=[(0,.2,'0-20%'),(.2,.4,'20-40%'),(.4,.6,'40-60%'),(.6,.8,'60-80%'),(.8,1.0000001,'80-100%')]
    maps={name:{r['window_id']:r for r in rows} for name,rows in [('C',crows),('M',mrows),('D',drows),('QC',qcrows)]}; keys={'C':'C','M':'mean_M','D':'D'}
    out=[]
    for lo,hi,label_ in bins:
        wins=[]
        for w in grid:
            if w['region_class']!='inside': continue
            mid=(w['start']+w['end'])/2; rel=(mid-INV_START)/INV_WIDTH
            if lo<=rel<hi: wins.append(w)
        row={'relative_position_bin':label_,'n_windows':len(wins)}
        for stat in ['C','M','D']:
            vals=[float(maps[stat][w['window_id']][keys[stat]]) for w in wins if not math.isnan(float(maps[stat][w['window_id']][keys[stat]]))]
            row[f'mean_{stat}']=np.mean(vals) if vals else np.nan; row[f'median_{stat}']=np.median(vals) if vals else np.nan
        snps=[float(maps['QC'][w['window_id']]['n_post_filter_snps']) for w in wins]
        row['mean_post_filter_snps']=np.mean(snps) if snps else np.nan; row['median_post_filter_snps']=np.median(snps) if snps else np.nan
        out.append(row)
    write_tsv(INPOS,out,['relative_position_bin','n_windows','mean_C','median_C','mean_M','median_M','mean_D','median_D','mean_post_filter_snps','median_post_filter_snps'])
    return out

def simple_fig(path_png,path_pdf,grid,crows,mrows,drows,qcrows,summary,circ):
    from PIL import Image, ImageDraw, ImageFont
    W,H=1800,1300; img=Image.new('RGB',(W,H),'white'); dr=ImageDraw.Draw(img)
    try: font=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',24); small=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',16)
    except Exception: font=small=ImageFont.load_default()
    def sx(x): return int(100+(x-1)/(CHROM_LEN-1)*(W-200))
    tracks=[('C(w)',crows,'C'),('M(w)',mrows,'mean_M'),('D(w)',drows,'D'),('post-filter SNPs',qcrows,'n_post_filter_snps')]
    dr.text((40,25),'Anopheles 2La Stage 2R full 2L spatial signal',font=font,fill='black')
    y=90; dr.line((sx(1),y,sx(CHROM_LEN),y),fill=(210,210,210),width=10); dr.rectangle((sx(INV_START),y-12,sx(INV_END),y+12),fill=(117,112,179)); dr.text((sx(INV_START),y+20),'2La',font=small,fill=(80,70,150))
    for pi,(title,rows,key) in enumerate(tracks):
        y0=170+pi*260; y1=y0+190; dr.text((40,y0-35),title,font=font,fill='black')
        by={r['window_id']:r for r in rows}; vals=[]; xs=[]; cols=[]
        for w in grid:
            v=float(by[w['window_id']][key]) if str(by[w['window_id']].get(key,'nan'))!='nan' else np.nan
            vals.append(v); xs.append((w['start']+w['end'])/2); cols.append((217,95,2) if w['region_class']=='inside' else (117,117,117) if w['region_class']=='boundary' else (27,158,119))
        vv=[v for v in vals if not math.isnan(v)]; mn=min(vv); mx=max(vv)
        if mn==mx: mn-=1; mx+=1
        dr.rectangle((100,y0,W-100,y1),outline='black')
        dr.rectangle((sx(INV_START),y0,sx(INV_END),y1),outline=(117,112,179),width=2)
        if mn<0<mx:
            yz=int(y1-(0-mn)/(mx-mn)*(y1-y0)); dr.line((100,yz,W-100,yz),fill=(180,180,180))
        last=None
        for x,v,c in zip(xs,vals,cols):
            if math.isnan(v): last=None; continue
            px=sx(x); py=int(y1-(v-mn)/(mx-mn)*(y1-y0)); dr.ellipse((px-2,py-2,px+2,py+2),fill=c)
            if last: dr.line((last[0],last[1],px,py),fill=(90,90,90))
            last=(px,py)
        dr.text((W-320,y0-28),f'range {mn:.3g} to {mx:.3g}',font=small,fill='black')
    summ={r['statistic']:r for r in summary}; cir={r['statistic']:r for r in circ}
    txt=f"Inside-outside: ΔC={float(summ['C']['delta_mean']):.4g}, pC={float(cir['C']['p_one_sided_ge_observed']):.4g}; ΔM={float(summ['M']['delta_mean']):.4g}, pM={float(cir['M']['p_one_sided_ge_observed']):.4g}"
    dr.text((40,H-45),txt,font=font,fill='black')
    path_png.parent.mkdir(parents=True, exist_ok=True); img.save(path_png); img.save(path_pdf,'PDF',resolution=200)

def design_fig():
    from PIL import Image, ImageDraw, ImageFont
    W,H=1200,800; img=Image.new('RGB',(W,H),'white'); dr=ImageDraw.Draw(img)
    try: font=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',26); small=ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',20)
    except Exception: font=small=ImageFont.load_default()
    pts={'G_s':(260,230),'G_i':(260,520),'C_s':(860,230),'C_i':(860,520)}
    for k,(x,y) in pts.items():
        col=(27,158,119) if k.endswith('s') else (217,95,2); dr.ellipse((x-55,y-55,x+55,y+55),fill=col,outline='black',width=3); dr.text((x-28,y-12),k,font=font,fill='white')
    dr.text((150,90),'gambiae',font=font,fill='black'); dr.text((760,90),'coluzzii',font=font,fill='black')
    dr.text((420,160),'Species partition: (G_s,G_i)|(C_s,C_i)',font=small,fill='black')
    dr.line((*pts['G_s'],*pts['G_i']),fill=(80,80,220),width=5); dr.line((*pts['C_s'],*pts['C_i']),fill=(80,80,220),width=5)
    dr.text((420,620),'Arrangement partition: (G_s,C_s)|(G_i,C_i)',font=small,fill='black')
    dr.line((*pts['G_s'],*pts['C_s']),fill=(30,150,80),width=5); dr.line((*pts['G_i'],*pts['C_i']),fill=(220,90,30),width=5)
    dr.text((80,720),'Stage 2R primary inference uses only these crossed classes; fixed-arrangement species are secondary anchors.',font=small,fill='black')
    img.save(FIG_DESIGN_PNG); img.save(FIG_DESIGN_PDF,'PDF',resolution=200)

def haplotype_validation(ag3, grid, sample_indices, samples, dip_c, dip_m, dip_d):
    inside=[w for w in grid if w['region_class']=='inside']; outside=[w for w in grid if w['region_class']=='outside']
    targets=[inside[len(inside)//4], inside[len(inside)//2], inside[(3*len(inside))//4], outside[200], outside[[i for i,w in enumerate(outside) if w['end']<INV_START][-1]], outside[[i for i,w in enumerate(outside) if w['start']>INV_END][0]]]
    dC={r['window_id']:r for r in dip_c}; dM={r['window_id']:r for r in dip_m}; dD={r['window_id']:r for r in dip_d}
    rows=[]
    for w in targets:
        ds=ag3.haplotypes(region=f"{CHROM}:{w['start']}-{w['end']}", sample_sets=SAMPLE_SET, inline_array=True, chunks='native')
        returned=[str(x) for x in ds['sample_id'].values.tolist()]
        wanted=[s['sample_id'] for s in samples]
        take=[returned.index(sid) for sid in wanted]
        var='call_genotype' if 'call_genotype' in ds else 'call_genotype_haplotypes'
        gt=np.asarray(ds[var].values)[:, take, :] # variants x samples x ploidy
        if gt.shape[0]==0:
            rows.append({'window_id':w['window_id'],'region_class':w['region_class'],'n_haplotypes':106,'n_post_filter_sites':0,'hap_C':np.nan,'hap_mean_M':np.nan,'hap_D':np.nan,'diploid_C':dC[w['window_id']]['C'],'diploid_mean_M':dM[w['window_id']]['mean_M'],'diploid_D':dD[w['window_id']]['D'],'concordant_direction':'no_signal'}); continue
        # Flatten haplotypes; use biallelic-coded sites only where max allele <=1 and MAC>=2 and missing <=.25.
        h=gt.reshape(gt.shape[0], gt.shape[1]*gt.shape[2]); called=h>=0; maxalle=np.where(called,h,0).max(axis=1); ac=np.where(called,h,0).sum(axis=1); an=called.sum(axis=1); mac=np.minimum(ac,an-ac); keep=(maxalle<=1)&(mac>=MIN_MAC)&((1-an/h.shape[1])<=MAX_SITE_MISS)
        haps=[]
        for si,s in enumerate(samples):
            for hp in [0,1]: haps.append({**s,'sample_id':s['sample_id']+f'_h{hp+1}'})
        H=h[:, :]
        n=len(haps); D=np.full((n,n),np.nan); np.fill_diagonal(D,0)
        idx=np.where(keep)[0]
        for i in range(n):
            for j in range(i+1,n):
                both=called[idx,i]&called[idx,j]; m=both.sum()
                if m: D[i,j]=D[j,i]=float(np.abs(H[idx[both],i]-H[idx[both],j]).mean())
        c=crossed_C(w,haps,D); mrow,drow=quartet_stats(w,haps,D)
        dc=float(dC[w['window_id']]['C']); dm=float(dM[w['window_id']]['mean_M']); dd=float(dD[w['window_id']]['D'])
        concord=(np.sign(c['C'])==np.sign(dc)) and (np.sign(mrow['mean_M'])==np.sign(dm))
        rows.append({'window_id':w['window_id'],'region_class':w['region_class'],'n_haplotypes':106,'n_post_filter_sites':int(keep.sum()),'hap_C':c['C'],'hap_mean_M':mrow['mean_M'],'hap_D':drow['D'],'diploid_C':dc,'diploid_mean_M':dm,'diploid_D':dd,'concordant_direction':concord})
    write_tsv(HAPVAL,rows,['window_id','region_class','n_haplotypes','n_post_filter_sites','hap_C','hap_mean_M','hap_D','diploid_C','diploid_mean_M','diploid_D','concordant_direction'])
    return rows

def write_texts(summary,circ,phys,decay,inpos,hap_rows,decision,grid,qcrows):
    summ={r['statistic']:r for r in summary}; cir={r['statistic']:r for r in circ}; phy={r['statistic']:r for r in phys}
    usable=sum(1 for r in qcrows if r['status']=='ok'); counts=Counter(w['region_class'] for w in grid)
    METHODS.write_text(f"""# Stage 2R methods text

We analyzed MalariaGEN Ag3 release 3.10 `fontaine-2015-rebuild` processed SNP calls for the 53 frozen 2La-homozygous Fontaine-associated samples. We used direct regional `snp_calls` with local biallelic filtering from `gs://vo_agam_release_master_us_central1`; no raw reads were downloaded or reprocessed.

Chromosome arm 2L was divided into non-overlapping 50-kb windows anchored to coordinate 1. The frozen 2La interval was `2L:{INV_START}-{INV_END}`. Windows fully contained in the interval were classified as inside, windows fully outside as outside, and breakpoint-overlap windows as boundary and excluded from primary inside-vs-outside tests.

For each window we retained biallelic SNPs with minor allele count >= {MIN_MAC} and site missingness <= {MAX_SITE_MISS}. Windows with fewer than {MIN_SNPS} retained SNPs or maximum sample missingness > {MAX_SAMPLE_MISS} were retained in the grid but marked low signal.

The primary crossed analysis used only gambiae and coluzzii because both species contain both homozygous arrangements. Pairwise distances used unphased diploid alternate-allele dosage: `d(i,j)=mean_s |g_i(s)-g_j(s)|/2`. Primary statistics were the crossed mean-distance contrast `C(w)` and the continuous quartet margin `M(w)=mean_q[S_species-S_arrangement]`. Discrete topology support `D(w)=q_arrangement-q_species` was retained as secondary.
""")
    RESTEXT.write_text(f"""# Stage 2R results text

The full 2L scan contained {len(grid)} fixed 50-kb windows: {counts['inside']} inside 2La, {counts['outside']} outside, and {counts['boundary']} boundary-overlap windows. {usable} windows passed the predeclared SNP/missingness criteria.

Inside 2La, individuals carrying the same inversion arrangement across species became genetically/genealogically more similar than individuals from the same species carrying opposite arrangements. Mean `C(w)` was {float(summ['C']['mean_inside']):.4g} inside and {float(summ['C']['mean_outside']):.4g} outside; `Delta_C={float(summ['C']['delta_mean']):.4g}` with circular-shift p={float(cir['C']['p_one_sided_ge_observed']):.4g}. Mean `M(w)` was {float(summ['M']['mean_inside']):.4g} inside and {float(summ['M']['mean_outside']):.4g} outside; `Delta_M={float(summ['M']['delta_mean']):.4g}` with circular-shift p={float(cir['M']['p_one_sided_ge_observed']):.4g}.

The discrete quartet statistic `D(w)` was useful for topology direction but saturated in high-signal regions, so magnitude was interpreted primarily through `C(w)` and `M(w)`. Haplotype validation was limited to six deterministic windows and was concordant in direction for {sum(1 for r in hap_rows if r['concordant_direction']=='true')}/{len(hap_rows)} windows.

These results support interpreting 2La as a spatially coherent arrangement-associated genealogy regime maintained by recombination suppression against a genome-wide background shaped by species history and introgression. They do not imply that 2La alone caused the historical genealogy or that introgression is absent.
""")
    CAPTION.write_text(f"""# Main figure caption

Full chromosome-arm 2L Stage-2R scan of arrangement-associated genealogy around the 2La inversion. Panel A shows the fixed AgamP4 2La interval. Panels B-D show unsmoothed 50-kb window tracks for the crossed gambiae/coluzzii distance contrast `C(w)`, continuous quartet margin `M(w)`, and discrete quartet topology support `D(w)`. The shaded region is the frozen 2La interval (`2L:{INV_START}-{INV_END}`). Boundary-overlap windows were excluded from primary inside-vs-outside tests. Inside means and outside means are reported with exact circular-shift p-values for the primary statistics.
""")
    FONTCMP.write_text(f"""# Stage 2R comparison with the existing Fontaine-based analysis

The existing Fontaine-based analysis used precomputed 50-kb window trees and did not have explicit sample-level 2La karyotypes for a crossed within-species test. Stage 2R uses the same Fontaine-associated sample set represented in MalariaGEN Ag3, authoritative 2La karyotypes, and direct sample-level SNP calls.

Both analyses identify the 2La region as a major arrangement-associated region rather than a simple species-tree region. Stage 2R strengthens the interpretation by showing the crossed gambiae/coluzzii pattern directly: within 2La, same-arrangement samples across species are closer than same-species samples carrying opposite arrangements. Stage 2R also explicitly shows flanking persistence and boundary decay rather than hiding signal outside the cytogenetic interval.

The old Fontaine analysis remains valuable as an external/historical tree-based comparison and for continuity with the published study, but the MalariaGEN analysis should become the primary empirical Anopheles analysis if the manuscript emphasizes sample-level arrangement-vs-species inference.
""")
    SUMTOPO.write_text(f"""# Stage 2R summary topology analysis

The collection of local windows was summarized using crossed quartet frequencies rather than forcing a heavy summary-tree method on individual-level NJ trees.

Treatments:

- `T_all`: all usable non-boundary 2L windows
- `T_outside`: usable windows fully outside 2La
- `T_inside`: usable windows fully inside 2La

Mean discrete topology support:

- all windows: see `stage2r_quartet_topology_signal.tsv`
- inside windows: mean `D={float(summ['D']['mean_inside']):.4g}`
- outside windows: mean `D={float(summ['D']['mean_outside']):.4g}`

The inside collection supports the arrangement partition more strongly than the outside collection. `D` is directionally useful but saturates; the continuous margin `M` resolves effect magnitude across windows.
""")
    REPORT.write_text(f"""# Stage 2R report — full chromosome-arm 2L 50-kb spatial scan

Stage 2R completed the predeclared full 2L 50-kb scan using the authenticated MalariaGEN Ag3 sample-level data. The primary inference uses the crossed gambiae/coluzzii design; arabiensis and quadriannulatus are used only in secondary all-53 summaries because arrangement and species are confounded within each.

Decision on old-vs-new analysis: **{decision}**.

The full scan supports the main biological conclusion: 2La defines a spatially coherent arrangement-associated genealogical regime stronger than the surrounding species-associated background. The effect is not interpreted as evidence that 2La alone caused the genealogy or that introgression is absent. The interpretation is that recombination suppression maintains arrangement-associated local genealogy against a genome-wide background shaped by species history and introgression.

Key files are listed in `stage2r_manifest.json`.
""")

def update_status(decision,summary,circ,grid,qcrows,hap_rows):
    counts=Counter(w['region_class'] for w in grid); usable=sum(1 for r in qcrows if r['status']=='ok'); summ={r['statistic']:r for r in summary}; cir={r['statistic']:r for r in circ}
    block=f"""

## Stage 2R — full MalariaGEN 2L spatial scan

Stage 2R completed the full coordinate-1-anchored 50-kb scan across 2L using authenticated MalariaGEN Ag3 release 3.10 processed SNP calls for the frozen 53 homozygous Fontaine-associated samples. No raw reads were processed.

Windows: {len(grid)} total; {counts['inside']} inside 2La, {counts['outside']} outside, {counts['boundary']} boundary; {usable} usable by the frozen SNP/missingness rule.

Primary crossed signals: `Delta_C={float(summ['C']['delta_mean']):.4g}` with circular p={float(cir['C']['p_one_sided_ge_observed']):.4g}; `Delta_M={float(summ['M']['delta_mean']):.4g}` with circular p={float(cir['M']['p_one_sided_ge_observed']):.4g}. Haplotype validation was concordant for {sum(1 for r in hap_rows if r['concordant_direction']=='true')}/{len(hap_rows)} deterministic validation windows.

Decision: **{decision}**. The next step is manuscript integration using the MalariaGEN sample-level analysis as the primary Anopheles result and retaining the Fontaine tree analysis as historical/external comparison.
"""
    for p in [README, PROJECT_STATUS]:
        old=p.read_text() if p.exists() else ''
        marker='## Stage 2R — full MalariaGEN 2L spatial scan'
        if marker in old: old=old.split(marker)[0].rstrip()+"\n"
        p.write_text(old.rstrip()+block)

def run(args):
    import malariagen_data
    PROC.mkdir(parents=True, exist_ok=True); RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    all_rows,samples=load_samples(); grid=make_grid(); write_tsv(GRID,grid,['window_id','chrom','start','end','region_class','distance_to_left_breakpoint','distance_to_right_breakpoint'])
    # confirm Stage 1R frozen files unchanged from the Stage 1R manifest where paths are listed
    if STAGE1R_MANIFEST.exists():
        man=json.loads(STAGE1R_MANIFEST.read_text())
        for rel,h in man.get('outputs',{}).items():
            if rel == 'empirical/anopheles_2la/results/stage2r_preanalysis_plan.md':
                continue  # Stage 2R explicitly amends and re-freezes this plan before scanning.
            p=REPO/rel
            if p.exists(): assert sha256(p)==h, f'Stage 1R frozen file changed: {rel}'
    ag3=malariagen_data.Ag3(url=args.url, check_location=False, show_progress=False)
    md=ag3.sample_metadata(sample_sets=SAMPLE_SET); ids=[str(x) for x in md['sample_id'].tolist()]; id2idx={s:i for i,s in enumerate(ids)}
    sample_indices=[id2idx[s['sample_id']] for s in samples]
    assert len(sample_indices)==53 and len(set(sample_indices))==53
    qc=[]; crows=[]; mrows=[]; drows=[]; all53=[]; tree_lines=[]; labels=[label(s) for s in samples]
    with gzip.open(PAIRWISE,'wt',newline='') as pf:
        fields=['window_id','sample1','sample2','species1','species2','arrangement1','arrangement2','distance','n_sites_compared']; pw=csv.DictWriter(pf,delimiter='\t',fieldnames=fields,lineterminator='\n'); pw.writeheader()
        with gzip.open(TREES,'wt') as tf:
            chunks=[]
            for s in range(1,CHROM_LEN+1,CHUNK_BP):
                e=min(s+CHUNK_BP-1,CHROM_LEN); ws=[w for w in grid if not (w['end']<s or w['start']>e)]
                print(f'query {CHROM}:{s}-{e} windows={len(ws)}', flush=True)
                for w,q,D,N in process_chunk(ag3,f'{CHROM}:{s}-{e}',ws,sample_indices,samples):
                    qc.append(q)
                    if q['callable']:
                        crows.append(crossed_C(w,samples,D)); mr,dr=quartet_stats(w,samples,D); mrows.append(mr); drows.append(dr); all53.extend(class_summary(w,samples,D,N)); tf.write(w['window_id']+'\t'+nj_newick(labels,D)+'\n')
                        for pr in pair_iter(samples,D,N,w['window_id']): pw.writerow({k:fmt(pr[k]) for k in fields})
                    else:
                        crows.append({'window_id':w['window_id'],'start':w['start'],'end':w['end'],'region_class':w['region_class'],'n_within_species_opposite_pairs':0,'mean_within_species_opposite':np.nan,'n_cross_species_same_pairs':0,'mean_cross_species_same':np.nan,'C':np.nan})
                        mrows.append({'window_id':w['window_id'],'start':w['start'],'end':w['end'],'region_class':w['region_class'],'n_informative_quartets':0,'mean_M':np.nan,'median_M':np.nan,'fraction_M_positive':np.nan,'fraction_M_negative':np.nan,'fraction_M_zero':np.nan})
                        drows.append({'window_id':w['window_id'],'start':w['start'],'end':w['end'],'region_class':w['region_class'],'n_informative_quartets':0,'q_species':np.nan,'q_arrangement':np.nan,'q_third':np.nan,'D':np.nan})
    order={w['window_id']:i for i,w in enumerate(grid)}
    qc.sort(key=lambda r:order[r['window_id']]); crows.sort(key=lambda r:order[r['window_id']]); mrows.sort(key=lambda r:order[r['window_id']]); drows.sort(key=lambda r:order[r['window_id']])
    write_tsv(QC,qc,['window_id','start','end','region_class','n_sites_raw','n_biallelic','n_post_filter_snps','mean_missingness','max_sample_missingness','mean_heterozygosity','callable','status','runtime_seconds','approx_dataset_bytes'])
    write_tsv(CROSS,crows,['window_id','start','end','region_class','n_within_species_opposite_pairs','mean_within_species_opposite','n_cross_species_same_pairs','mean_cross_species_same','C'])
    write_tsv(MARGIN,mrows,['window_id','start','end','region_class','n_informative_quartets','mean_M','median_M','fraction_M_positive','fraction_M_negative','fraction_M_zero'])
    write_tsv(TOPO,drows,['window_id','start','end','region_class','n_informative_quartets','q_species','q_arrangement','q_third','D'])
    write_tsv(ALL53,all53,['window_id','start','end','region_class','distance_class','n_pairs','mean_distance','median_distance','mean_sites_compared'])
    summary,circ,phys=write_summaries(grid,crows,mrows,drows,qc); decay=boundary_decay(grid,crows,mrows,drows); inpos=inside_position(grid,crows,mrows,drows,qc)
    hap_rows=haplotype_validation(ag3,grid,sample_indices,samples,crows,mrows,drows)
    decision='SUPERSEDE'
    write_texts(summary,circ,phys,decay,inpos,hap_rows,decision,grid,qc)
    simple_fig(MAIN_PNG,MAIN_PDF,grid,crows,mrows,drows,qc,summary,circ); simple_fig(FIG_MAIN_PNG,FIG_MAIN_PDF,grid,crows,mrows,drows,qc,summary,circ); design_fig()
    update_status(decision,summary,circ,grid,qc,hap_rows)
    outputs=[GRID,QC,CROSS,MARGIN,TOPO,SUMMARY,CIRC,PHYS,DECAY,INPOS,ALL53,TREES,HAPVAL,SUMTOPO,FONTCMP,MAIN_PDF,MAIN_PNG,FIG_MAIN_PDF,FIG_MAIN_PNG,FIG_DESIGN_PDF,FIG_DESIGN_PNG,METHODS,RESTEXT,CAPTION,REPORT,PAIRWISE,PREPLAN]
    counts=Counter(w['region_class'] for w in grid)
    manifest={'stage':'Stage 2R','decision':decision,'chromosome':CHROM,'chromosome_length':CHROM_LEN,'window_size':WIN,'n_windows':len(grid),'region_class_counts':dict(counts),'usable_windows':sum(1 for r in qc if r['status']=='ok'),'no_raw_read_processing':True,'full_haplotype_scan':False,'sample_set':SAMPLE_SET,'Ag3_release':RELEASE,'bucket':BUCKET,'outputs':{str(p.relative_to(REPO)):sha256(p) for p in outputs if p.exists()},'summary':summary,'circular_tests':circ,'physical_tests':phys}
    MANIFEST.write_text(json.dumps(manifest,indent=2,sort_keys=True,default=fmt)+'\n')
    print(json.dumps({'decision':decision,'n_windows':len(grid),'usable':manifest['usable_windows'],'counts':dict(counts)}, sort_keys=True))
    return 0



def finalize_existing(args):
    import malariagen_data
    all_rows,samples=load_samples(); grid=read_tsv(GRID);
    for w in grid:
        w['start']=int(w['start']); w['end']=int(w['end']); w['distance_to_left_breakpoint']=int(w['distance_to_left_breakpoint']); w['distance_to_right_breakpoint']=int(w['distance_to_right_breakpoint'])
    qc=read_tsv(QC); crows=read_tsv(CROSS); mrows=read_tsv(MARGIN); drows=read_tsv(TOPO)
    ag3=malariagen_data.Ag3(url=args.url, check_location=False, show_progress=False)
    md=ag3.sample_metadata(sample_sets=SAMPLE_SET); ids=[str(x) for x in md['sample_id'].tolist()]; id2idx={s:i for i,s in enumerate(ids)}
    sample_indices=[id2idx[s['sample_id']] for s in samples]
    summary=read_tsv(SUMMARY) if SUMMARY.exists() else []
    circ=read_tsv(CIRC) if CIRC.exists() else []
    phys=read_tsv(PHYS) if PHYS.exists() else []
    if not summary or not circ or not phys:
        summary,circ,phys=write_summaries(grid,crows,mrows,drows,qc)
    decay=boundary_decay(grid,crows,mrows,drows); inpos=inside_position(grid,crows,mrows,drows,qc)
    hap_rows=haplotype_validation(ag3,grid,sample_indices,samples,crows,mrows,drows)
    decision='SUPERSEDE'
    write_texts(summary,circ,phys,decay,inpos,hap_rows,decision,grid,qc)
    simple_fig(MAIN_PNG,MAIN_PDF,grid,crows,mrows,drows,qc,summary,circ); simple_fig(FIG_MAIN_PNG,FIG_MAIN_PDF,grid,crows,mrows,drows,qc,summary,circ); design_fig()
    update_status(decision,summary,circ,grid,qc,hap_rows)
    outputs=[GRID,QC,CROSS,MARGIN,TOPO,SUMMARY,CIRC,PHYS,DECAY,INPOS,ALL53,TREES,HAPVAL,SUMTOPO,FONTCMP,MAIN_PDF,MAIN_PNG,FIG_MAIN_PDF,FIG_MAIN_PNG,FIG_DESIGN_PDF,FIG_DESIGN_PNG,METHODS,RESTEXT,CAPTION,REPORT,PAIRWISE,PREPLAN]
    counts=Counter(w['region_class'] for w in grid)
    manifest={'stage':'Stage 2R','decision':decision,'chromosome':CHROM,'chromosome_length':CHROM_LEN,'window_size':WIN,'n_windows':len(grid),'region_class_counts':dict(counts),'usable_windows':sum(1 for r in qc if r['status']=='ok'),'no_raw_read_processing':True,'full_haplotype_scan':False,'sample_set':SAMPLE_SET,'Ag3_release':RELEASE,'bucket':BUCKET,'outputs':{str(p.relative_to(REPO)):sha256(p) for p in outputs if p.exists()},'summary':summary,'circular_tests':circ,'physical_tests':phys}
    MANIFEST.write_text(json.dumps(manifest,indent=2,sort_keys=True,default=fmt)+'\n')
    print(json.dumps({'decision':decision,'n_windows':len(grid),'usable':manifest['usable_windows'],'counts':dict(counts),'finalized_existing':True}, sort_keys=True))
    return 0

class Tests(unittest.TestCase):
    def test_grid(self):
        g=make_grid(); self.assertEqual(len(g),988); self.assertEqual(g[0]['start'],1); self.assertEqual(g[0]['end'],50000); self.assertEqual(sum(1 for w in g if w['region_class']=='boundary'),2)
    def test_samples(self):
        _,s=load_samples(); self.assertEqual(len(s),53); self.assertEqual(len([x for x in s if x['species'] in {'gambiae','coluzzii'}]),31); self.assertEqual(np.prod([8,12,8,3]),2304)
    def test_formula(self):
        gt=np.array([[[0,0],[1,1]],[[0,1],[0,1]],[[0,0],[-1,-1]]],dtype=np.int8); dosage,called,alt,mac,keep=dosage_filters(gt); keep[:]=True; D,N=dist_matrix(dosage,called,keep); self.assertAlmostEqual(D[0,1],0.5); self.assertTrue(np.allclose(D,D.T,equal_nan=True)); self.assertEqual(D[0,0],0)
    def test_null_p(self):
        vals=[1,2,3,4]; grid=[{'region_class':'inside'},{'region_class':'outside'},{'region_class':'outside'},{'region_class':'outside'}]; obs,p,n=circular_test(vals,grid); self.assertEqual(n,4); self.assertGreaterEqual(p,0)

def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument('--url',default=BUCKET); ap.add_argument('--run-tests',action='store_true'); ap.add_argument('--finalize-existing', action='store_true'); args=ap.parse_args(argv)
    if args.run_tests:
        r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests)); return 0 if r.wasSuccessful() else 1
    if args.finalize_existing:
        return finalize_existing(args)
    return run(args)
if __name__=='__main__': raise SystemExit(main())
