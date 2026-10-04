#!/usr/bin/env python3
"""Exhaustive named/collapsed/coordinate audit, before any inference or freeze."""
import csv
import gzip
import json
import lzma
import re
from collections import Counter, defaultdict
from itertools import zip_longest
from stage4d_core import *

def main():
    stage4a = load_stage('04a_structural_topology_support.py')
    masks = stage4a.build_masks(stage4a.read_tsv(stage4a.BREAKPOINTS))
    meta_path = ROOT / 'data/neoaves_chr4/raw/genetreesupport/63K_trees.names_header.txt.xz'
    with lzma.open(meta_path, 'rt') as f:
        meta = list(csv.DictReader(f, delimiter=' ', skipinitialspace=True))
    by_key = {}
    for r in meta:
        key = (r['Chromosome'], int(r['ws']), int(r['we']), int(r['wstart']))
        if key in by_key: raise ValueError('duplicate coordinate metadata')
        by_key[key] = r
    master = {}
    with (EXTERNAL / 'master_table_gene_trees.txt').open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['part_of_63K'] == 'Y':
                if r['locus'] in master: raise ValueError('duplicate master locus')
                master[r['locus']] = r
    import dendropy
    pub = dendropy.Tree.get(path=str(EXTERNAL / '63K.tre'), schema='newick', preserve_underscores=True)
    taxa = sorted(n.taxon.label for n in pub.leaf_node_iter())
    bits = {taxon: 1 << i for i, taxon in enumerate(taxa)}
    named = EXTERNAL / '63430.named.gene.trees.gz'
    collapsed = EXTERNAL / '63430.aLRT-0.95-collapsed.gene.trees.gz'
    rows, ids, mapped = [], set(), set()
    completeness = Counter()
    audit = Counter()
    max_depth_error = 0
    signatures = Counter()
    DATA.mkdir(parents=True, exist_ok=True)
    original_input = DATA / 'FULL363_ORIGINAL_REPRO.tre'
    mismatches=[]
    # The collapsed file has no IDs. Test the documented row-order claim
    # exhaustively, but never assign IDs by row after the first disagreement.
    with gzip.open(named, 'rt') as nf, gzip.open(collapsed, 'rt') as cf, original_input.open('w') as archived_out, (DATA / 'named_recollapsed.tre').open('w') as out:
        for index, (nl,cl) in enumerate(zip_longest(nf,cf), 1):
            if nl is None or cl is None: raise ValueError('named/collapsed counts differ')
            archived_out.write(cl)
            locus, nwk = nl.strip().split(None, 1)
            if locus in ids: raise ValueError('duplicate locus ID')
            ids.add(locus)
            m = re.fullmatch(r'(.+)_(\d+)_(\d+)_1k_start(\d+)_fasta.treefile', locus)
            if not m: raise ValueError(f'unknown locus schema: {locus}')
            chrom, ws, we, offset = m.groups()
            ws, we, offset = int(ws), int(we), int(offset)
            r = by_key[(chrom, ws, we, offset)]
            mapped.add(r['Gene'])
            mr = master[locus.removesuffix('_fasta.treefile')]
            if mr['start'] != 'NA':
                assert (mr['chromosome'], int(float(mr['start'])), int(float(mr['end']))) == (chrom, ws, we)
                audit['master_coordinates_verified'] += 1
            else:
                audit['master_coordinates_missing'] += 1
            parsed = parse_newick(nwk)
            a = tree_signature(parsed, bits, .95)
            b = tree_signature(parse_newick(cl.strip()), bits)
            row_match = a[0] == b[0] and a[1] == b[1]
            if row_match:
                audit['same_row_taxa_and_collapsed_splits'] += 1
                error=max(abs(a[3][t]-b[3][t]) for t in a[3])
                max_depth_error=max(max_depth_error,error)
                audit['same_row_path_lengths_within_1e-4'] += error <= .0001
            else:
                mismatches.append(dict(tree_index=index,locus_id=locus,
                    named_n_taxa=len(a[2]),collapsed_same_row_n_taxa=len(b[2]),
                    taxon_set_matches=a[0]==b[0],split_set_matches=a[1]==b[1]))
            # Full precision named leaf paths distinguish duplicate topology trees.
            sig = (a[0], tuple(sorted(a[1])), tuple((t, round(a[3][t], 7)) for t in a[2]))
            signatures[sha256_bytes(repr(sig).encode())] += 1
            if len(a[2]) != int(mr['nTaxa']): raise ValueError('taxon count vs metadata')
            completeness[len(a[2])] += 1
            midpoint = (ws + we) / 2
            flags = {s: chrom == 'chr4' and stage4a.in_intervals(midpoint, masks[s]) for s in masks}
            rows.append(dict(locus_id=locus, tree_index=index, collapsed_tree_index='', original_gene_id=r['Gene'], chromosome=chrom,
                             start=ws, end=we, midpoint=midpoint, locus_start_within_window=offset,
                             n_taxa=len(a[2]), original_tree_file=str(named.relative_to(ROOT)),
                             collapsed_tree_file=str(collapsed.relative_to(ROOT)), is_chr4=is_chr4(chrom),
                             published_outlier_region=published_outlier(chrom, ws),
                             frozen_structural_stringent=flags['stringent'], frozen_structural_primary=flags['primary'],
                             frozen_structural_inclusive=flags['inclusive'],
                             structural_mask_class=';'.join(s for s in masks if flags[s]) or ('unmapped_scaffold' if is_chr4(chrom) and chrom != 'chr4' else 'background')))
            rebuilt=collapse_named(parsed)
            out.write(emit_newick(rebuilt)+'\n')
            if index % 5000 == 0: print(f'audited {index}', flush=True)
    assert len(rows) == len(meta) == len(master) == len(mapped) == 63430
    frozen = read_tsv(ROOT / 'data/neoaves_chr4/processed/locus_table_chr4.tsv')
    lookup = {r['original_gene_id']: r for r in rows}
    # Compare the shared window-coordinate convention; do not substitute the
    # selected 1kb alignment's start for the frozen 10kb-window midpoint.
    for old in frozen:
        new = lookup[old['Gene']]
        assert (new['chromosome'], new['start'], new['end'], new['midpoint']) == (old['Chromosome'], int(old['start']), int(old['end']), float(old['midpoint']))
    write_tsv(RESULTS / 'stage4d_locus_manifest.tsv', rows)
    write_tsv(RESULTS / 'stage4d/named_collapsed_mismatches.tsv', mismatches,
              ['tree_index','locus_id','named_n_taxa','collapsed_same_row_n_taxa','taxon_set_matches','split_set_matches'])
    summary = dict(n_loci=len(rows), n_taxa_union=len(taxa), checks=dict(audit),
                   max_root_to_tip_length_error=max_depth_error, duplicate_full_fingerprints=sum(v-1 for v in signatures.values()),
                   completeness=dict(sorted(completeness.items())), coordinate_mapping_failures=0,
                   chr4_sequence_counts=dict(Counter(r['chromosome'] for r in rows if r['is_chr4'])),
                   n_nonchr4=sum(not r['is_chr4'] for r in rows), n_pnas=sum(r['published_outlier_region'] for r in rows),
                   n_struct={s:sum(r['frozen_structural_'+s] for r in rows) for s in masks},
                   masks=masks, metadata_sha256=sha256(meta_path),
                   collapsed_unmatched_named_loci=len(mismatches),
                   mapping_solution='Exhaustive same-row comparison rejects the documented ordering where taxa/retained splits differ. No collapsed IDs were assigned by row. Fallback: all named trees re-collapsed at support<0.95 while locus IDs remain in the manifest; that reconstructed input must reproduce the published topology before use.',
                   coordinates='Existing Stage1-4C 10kb window start/end and midpoint; offset retained separately. Structural masks apply only to anchored chr4, not unplaced scaffold coordinates.')
    (RESULTS / 'stage4d/input_audit.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k != 'completeness'}, indent=2))

def sha256_bytes(b):
    import hashlib
    return hashlib.sha256(b).hexdigest()

if __name__ == '__main__': main()
