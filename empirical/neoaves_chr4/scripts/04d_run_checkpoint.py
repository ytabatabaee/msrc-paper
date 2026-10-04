#!/usr/bin/env python3
"""Run a reproduction checkpoint, preserving each attempt and its provenance.

No novel treatments are exposed here. A successful process is not itself a
scientific checkpoint pass: topology and focal-relationship validation follow.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess

from stage4d_core import (DATA, RESULTS, ROOT, sha256, rf_comparison,
                          focal_split_states, frozen_focal_groups)

SOFTWARE = DATA / 'software/ASTRAL-887669674518f38c13c108d006f936cd54d77467'
INPUTS = {
    'FULL363_ORIGINAL_REPRO': 'FULL363_ORIGINAL_REPRO.tre',
    'FULL363_RECOLLAPSED_VALIDATION': 'named_recollapsed.tre',
    'J48_ALL': 'J48_ALL.tre',
    'J48_PNAS_OUTLIER_EXCLUDED': 'J48_PNAS_OUTLIER_EXCLUDED.tre',
}


def command(java, software, source, output, heap, threads):
    return [java, f'-Xmx{heap}', f'-Djava.library.path={software / "lib"}',
            '-cp', os.pathsep.join((str(software / 'main'), str(software / 'lib/*'))),
            'phylonet.coalescent.CommandLine', '-C', '-T', str(threads),
            '-s', '692', '-t', '3', '-i', str(source), '-o', str(output)]


def completed(returncode, tree, log):
    return (returncode == 0 and tree.is_file() and tree.stat().st_size > 0
            and tree.read_text().strip().endswith(';')
            and 'ASTRAL finished in ' in log.read_text())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkpoint', choices=INPUTS)
    p.add_argument('--java', required=True, help='Explicit Java executable')
    p.add_argument('--heap', default='8g')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--software', type=Path, default=SOFTWARE)
    args = p.parse_args()
    if args.threads < 1:
        p.error('--threads must be positive')
    # Preserve the original inventory and validate before expensive computation.
    for row in json.loads((RESULTS / 'stage4d/frozen_prior_inventory.json').read_text()):
        if sha256(ROOT / row['path']) != row['sha256']:
            raise RuntimeError(f"Frozen checksum mismatch: {row['path']}")
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    attempt = RESULTS / 'stage4d/attempts' / f'{args.checkpoint}_{stamp}'
    attempt.mkdir(parents=True, exist_ok=False)
    source = DATA / INPUTS[args.checkpoint]
    tree, log = attempt / 'tree.nwk', attempt / 'stderr.log'
    cmd = command(args.java, args.software.resolve(), source, tree, args.heap, args.threads)
    software_hashes = {
        str(x.relative_to(args.software)): sha256(x)
        for x in sorted(args.software.rglob('*'))
        if x.is_file() and x.suffix in {'.class', '.jar', '.dylib', '.so', '.dll'}
    }
    record = dict(checkpoint=args.checkpoint, started_utc=stamp, command=cmd,
                  input=str(source.relative_to(ROOT)), input_sha256=sha256(source),
                  java_version=subprocess.run([args.java, '-version'], capture_output=True,
                                              text=True, check=True).stderr,
                  software_sha256=software_hashes, status='running',
                  scientific_gate_passed=False)
    provenance = attempt / 'provenance.json'
    provenance.write_text(json.dumps(record, indent=2) + '\n')
    print(attempt, flush=True)
    with (attempt / 'stdout.log').open('w') as stdout, log.open('w') as stderr:
        result = subprocess.run(cmd, stdout=stdout, stderr=stderr)
    record.update(returncode=result.returncode,
                  finished_utc=datetime.now(timezone.utc).isoformat(),
                  status='completed' if completed(result.returncode, tree, log) else 'failed',
                  stderr_sha256=sha256(log))
    if record['status'] == 'completed':
        record['tree_sha256'] = sha256(tree)
        record['focal_relationships'] = focal_split_states(tree.read_text(), frozen_focal_groups())
        if args.checkpoint.startswith('FULL363'):
            reference = ROOT / 'empirical/neoaves_chr4/external/stiller2024/63K.tre'
            record['published_tree_comparison'] = rf_comparison(tree.read_text(), reference.read_text())
            record['reference_sha256'] = sha256(reference)
        record['validation_note'] = 'Requires focal extraction and scientific gate review; not auto-promoted.'
        canonical = RESULTS / 'stage4d/trees' / f'{args.checkpoint}.nwk'
        canonical.parent.mkdir(parents=True, exist_ok=True)
        if not canonical.exists() or canonical.stat().st_size == 0:
            shutil.copyfile(tree, canonical)
            record['saved_tree'] = str(canonical.relative_to(ROOT))
        elif sha256(canonical) != sha256(tree):
            record['validation_note'] += ' Existing nonempty tree differs; preserved both attempts.'
    provenance.write_text(json.dumps(record, indent=2) + '\n')
    raise SystemExit(0 if record['status'] == 'completed' else 1)


if __name__ == '__main__':
    main()
