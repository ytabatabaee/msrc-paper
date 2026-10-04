#!/usr/bin/env python3
"""Build only the pre-freeze reproduction-checkpoint inputs."""
import json
from stage4d_core import *

def main():
    rows=read_tsv(RESULTS/'stage4d_locus_manifest.tsv')
    taxa={r['stiller2024_taxon'] for r in read_tsv(DATA/'jarvis48_taxa.tsv')}
    assert len(taxa)==48
    source=DATA/'named_recollapsed.tre'
    outputs={'J48_ALL':DATA/'J48_ALL.tre','J48_PNAS_OUTLIER_EXCLUDED':DATA/'J48_PNAS_OUTLIER_EXCLUDED.tre'}
    handles={k:p.open('w') for k,p in outputs.items()}; counts={k:0 for k in outputs}; observed=set()
    try:
        with source.open() as f:
            for row,line in zip(rows,f,strict=True):
                tree=prune_tree(parse_newick(line),taxa)
                if tree is None: continue
                leaves=set(tree_signature(tree,{t:1<<i for i,t in enumerate(sorted(taxa))})[2])
                observed |= leaves
                if len(leaves)>=4:
                    handles['J48_ALL'].write(emit_newick(tree)+'\n');counts['J48_ALL']+=1
                    if row['published_outlier_region']!='True':
                        handles['J48_PNAS_OUTLIER_EXCLUDED'].write(emit_newick(tree)+'\n');counts['J48_PNAS_OUTLIER_EXCLUDED']+=1
    finally:
        for h in handles.values():h.close()
    assert observed==taxa
    assert counts=={'J48_ALL':62945,'J48_PNAS_OUTLIER_EXCLUDED':61524},counts
    report={'jarvis_taxa':sorted(taxa),'counts':counts,'trees_with_fewer_than_four_jarvis_taxa_omitted':63430-counts['J48_ALL'],
            'usable_published_outlier_trees_omitted':counts['J48_ALL']-counts['J48_PNAS_OUTLIER_EXCLUDED'],
            'source':str(source.relative_to(ROOT)),
            'source_sha256':sha256(source),'taxon_mapping':str((DATA/'jarvis48_taxa.tsv').relative_to(ROOT))}
    (RESULTS/'stage4d/j48_input_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(counts))

if __name__=='__main__':main()
