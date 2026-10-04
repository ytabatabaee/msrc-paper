#!/usr/bin/env python3
"""Reproduce published Figure-4 summaries from authors' per-locus scores only.

This validates aggregation, not de-novo quartet scoring or exact taxon membership.
No novel MSRC treatment is evaluated here; those are gated on all checkpoints.
"""
import csv
import json
import lzma
import math
import re
from collections import defaultdict
import numpy as np
from stage4d_core import *

SOURCE = BASE / 'external/chr4avian'
DOVES = {'Alopecoenas_beccarii','Caloenas_nicobarica','Columba_livia','Columbina_picui','Patagioenas_fasciata'}
CUCKOOS = {'Centropus_bengalensis','Centropus_unirufus','Ceuthmochares_aereus','Crotophaga_sulcirostris',
           'Cuculus_canorus','Geococcyx_californianus','Piaya_cayana'}
ALIASES = {
 'Alopecoenas':'Alopecoenas_beccarii','Caloenas':'Caloenas_nicobarica','Columba':'Columba_livia',
 'Columbina':'Columbina_picui','Patagioenas':'Patagioenas_fasciata','Centropus_ben':'Centropus_bengalensis',
 'Centropus_uni':'Centropus_unirufus','Ceuthmochare':'Ceuthmochares_aereus','Crotophaga':'Crotophaga_sulcirostris',
 'Cuculus':'Cuculus_canorus','Geococcyx':'Geococcyx_californianus','Piaya_cayana':'Piaya_cayana',
 'Mesitornis_unicolor':'Mesitornis_unicolor','Podiceps_crus':'Podiceps_cristatus',
 'Podilymbus_podiceps':'Podilymbus_podiceps','Pterocles_burchelli':'Pterocles_burchelli',
 'Pterocles_gutturalis':'Pterocles_gutturalis','Syrrhaptes_paradoxus':'Syrrhaptes_paradoxus'}

def removed_taxa(name):
    if name == 'full-tree': return set()
    if name == 'all-doves': return DOVES
    if name == 'allcuckoo': return CUCKOOS
    if name == 'all-doves_Cuculus': return DOVES|{'Cuculus_canorus'}
    if name == 'alldoves-allcucko': return DOVES|CUCKOOS
    result=set(); remaining=name
    for alias in sorted(ALIASES,key=len,reverse=True):
        if alias in remaining:
            result.add(ALIASES[alias]); remaining=remaining.replace(alias,'')
    if remaining.strip('_-'): raise ValueError(f'unparsed removal condition {name}: {remaining}')
    return result

def condition_label(s):
    return re.sub(r'_([A-Z])', r' \1', s).replace('minus-', '')

def summarize(values):
    v = np.asarray(values, dtype=float)
    return dict(n_loci=len(v), mean=float(np.mean(v)), median=float(np.median(v)),
                se=float(np.std(v, ddof=1)/math.sqrt(len(v))),
                abs_mean_minus_median=float(abs(np.mean(v)-np.median(v))))

def main():
    rc = read_tsv(SOURCE / 'removed-count.tsv')
    # R factor/reorder sorts by the score; ties retain factor label order.
    rc.sort(key=lambda r: (100*int(r['Count'])+10*int(r['Dove'] or 0)+int(r['Cuckoos'] or 0), condition_label(r['Tree']).casefold()))
    with lzma.open(ROOT / 'data/neoaves_chr4/raw/genetreesupport/63K_trees.names_header.txt.xz', 'rt') as f:
        meta = list(csv.DictReader(f, delimiter=' ', skipinitialspace=True))
    outliers = {int(r['Gene']) for r in meta if published_outlier(r['Chromosome'], int(r['ws']))}
    values = defaultdict(list)
    genes = defaultdict(set)
    empty = defaultdict(int)
    with lzma.open(SOURCE / 'all.stat.xz', 'rt') as f:
        for line in f:
            name, gene, main, alt, diff = line.split()
            name = name.removeprefix('minus-')
            gene = int(gene); main, alt, diff = float(main), float(alt), float(diff)
            if abs(diff-(main-alt)) > 1e-6: raise ValueError('score sign inconsistency')
            if gene in genes[name]: raise ValueError('duplicate condition/locus')
            genes[name].add(gene)
            empty[name] += main + alt == 0
            values[(name, 'T0')].append(-diff)
            if gene not in outliers: values[(name, 'T_PNAS')].append(-diff)
    assert set(genes) == {r['Tree'] for r in rc}, (set(genes) - {r['Tree'] for r in rc}, {r['Tree'] for r in rc} - set(genes))
    rows=[]
    for order, r in enumerate(rc):
        for treatment in ('T0', 'T_PNAS'):
            rows.append(dict(condition=r['Tree'], label=condition_label(r['Tree']), axis_order=order,
                             n_removed=int(r['Count']), treatment=treatment,
                             **summarize(values[(r['Tree'], treatment)]),
                             source='authors all.stat.xz; aggregation reproduction only',
                             se_interpretation='original descriptive locus SE; not independence-valid MSRC uncertainty'))
    write_tsv(RESULTS / 'stage4d_figure4_reproduction.tsv', rows)
    grid=[]
    for i,r in enumerate(rc):
        removed=removed_taxa(r['Tree'])
        if len(removed) != int(r['Count']): raise ValueError((r['Tree'],len(removed),r['Count']))
        grid.append(dict(axis_order=i, **r, removed_taxa=';'.join(sorted(removed)),
                         taxon_membership_status='exact removal set deterministically expanded from authors condition label and verified against Count/Dove/Cuckoos columns'))
    write_tsv(RESULTS / 'stage4d/figure4_original_grid.tsv', grid)
    report = dict(n_conditions=len(rc), n_pnas_loci=len(outliers),
                  counts={k:len(v) for k,v in genes.items()}, zero_score_rows=dict(empty),
                  statistic='-(main-alt) = alt-main, as q$diff=-q$diff in delta_quartet.R',
                  standard_error='sample SD / sqrt(n); ggplot2 mean_se default mult=1',
                  taxon_membership_verified=True, de_novo_scores_verified=False,
                  axis_order='100*Count+10*Dove+Cuckoos, with original display-label factor tie ordering',
                  novel_treatments_run=False)
    (RESULTS / 'stage4d/figure4_reproduction_audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__': main()
