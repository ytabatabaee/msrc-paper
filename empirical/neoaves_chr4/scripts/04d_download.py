#!/usr/bin/env python3
"""Resolve selected Stiller resources from the saved index; resume, verify, manifest."""
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DEST = ROOT / 'empirical/neoaves_chr4/external/stiller2024'
BASE = 'https://erda.ku.dk/archives/341f72708302f1d0c461ad616e783b86/'
WANTED = [
 '01_alignments_and_gene_trees/intergenic_regions/63430.named.gene.trees.gz',
 '01_alignments_and_gene_trees/intergenic_regions/63430.aLRT-0.95-collapsed.gene.trees.gz',
 '03_species_trees/63K.tre',
 '02_characteristics_of_alignments_and_gene_trees/master_table_gene_trees.txt',
 '05_subsetting_experiments/chromosome/by_chromosome.tar.gz',
 '03_species_trees/ASTRAL_log_files/astral-63430_filter_mRNA_random_region.gene.trees.aLRT0.05.log',
]

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def main():
    index = json.loads((DEST / 'published-files.json').read_text())
    manifest = DEST / 'download_manifest.tsv'
    prior = {}
    if manifest.exists():
        with manifest.open() as f: prior = {r['filename']: r for r in csv.DictReader(f, delimiter='\t')}
    selected = [('published-archive.html', None), ('published-files.json', None)]
    for suffix in ['README.txt'] + WANTED:
        matches = [r for r in index if r['name'] == 'B10K/data_upload/' + suffix]
        if len(matches) != 1: raise ValueError(f'Archive resolution failed: {suffix}')
        selected.append((matches[0]['name'], matches[0]['size']))
    for archive_path, size in selected:
        p = DEST / Path(archive_path).name
        old = prior.get(p.name)
        if old and p.exists():
            if sha(p) != old['sha256']: raise ValueError(f'Immutable input changed: {p}')
            continue
        if not p.exists() or (size is not None and p.stat().st_size != size):
            subprocess.run(['curl', '-fL', '-C', '-', '--retry', '5', '--connect-timeout', '30', BASE + archive_path, '-o', str(p)], check=True)
        if size is not None and p.stat().st_size != size: raise ValueError(f'Partial download: {p}')
        integrity = 'size_verified' if size is not None else 'download_complete'
        if p.suffix == '.gz':
            subprocess.run(['gzip', '-t', str(p)], check=True)
            integrity += ';gzip_pass'
        if p.name.endswith('.tar.gz'):
            subprocess.run(['tar', '-tzf', str(p)], check=True, stdout=subprocess.DEVNULL)
            integrity += ';tar_pass'
        prior[p.name] = dict(filename=p.name, source_url=BASE + archive_path, archive_path=archive_path,
                            download_timestamp=datetime.now(timezone.utc).isoformat(), byte_size=p.stat().st_size,
                            sha256=sha(p), integrity=integrity)
        with manifest.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(prior[p.name]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(prior.values())
        p.chmod(0o444)
        print(p.name, integrity, flush=True)
    target = ROOT / 'empirical/neoaves_chr4/results/stage4d_download_manifest.tsv'
    target.write_bytes(manifest.read_bytes())

if __name__ == '__main__': main()
