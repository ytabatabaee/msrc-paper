#!/usr/bin/env python3
"""Supplementary recombination-suppression validation for Anopheles 2La.

This is not a new empirical stage. It adds an independent published crossing-data
recombination track alongside frozen Stage-2R C/M/D statistics.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, os, tempfile, textwrap, unittest
from pathlib import Path

import numpy as np

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'msrc-paper-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / 'data' / 'anopheles_2la'
RAW = DATA / 'raw' / 'recombination'
PROC = DATA / 'processed'
BASE = REPO / 'empirical' / 'anopheles_2la'
RESULTS = BASE / 'results'
FIGURES = BASE / 'figures'
README = BASE / 'README.md'

INV_START = 20_524_058
INV_END = 42_165_532
CHROM = '2L'
ASSEMBLY = 'AgamP4-style coordinates used by Stage 2R-4R'
SOURCE_CITATION = 'Stump AD, Pombi M, Goeddel L, Ribeiro JMC, Wilder JA, della Torre A, Besansky NJ. 2007. Genetic exchange in 2La inversion heterokaryotypes of Anopheles gambiae. Insect Molecular Biology 16(6):703-709.'
SOURCE_DOI = '10.1111/j.1365-2583.2007.00764.x'
SOURCE_URLS = [
    'https://doi.org/10.1111/j.1365-2583.2007.00764.x',
    'https://pubmed.ncbi.nlm.nih.gov/18092999/',
    'https://iris.uniroma1.it/handle/11573/363495',
    'https://experts.nau.edu/en/publications/genetic-exchange-in-2la-inversion-heterokaryotypes-of-anopheles-g/',
]

GRID = PROC / 'stage2r_window_grid.tsv'
CROSS = RESULTS / 'stage2r_crossed_distance_signal.tsv'
MARGIN = RESULTS / 'stage2r_quartet_margin_signal.tsv'
TOPO = RESULTS / 'stage2r_quartet_topology_signal.tsv'
STAGE2_MANIFEST = RESULTS / 'stage2r_manifest.json'
STAGE3_MANIFEST = RESULTS / 'stage3r_manifest.json'
STAGE4_MANIFEST = RESULTS / 'stage4r_manifest.json'

SOURCE_AUDIT = RESULTS / 'recombination_source_audit.md'
COORD_AUDIT = RESULTS / 'recombination_coordinate_audit.md'
SUMMARY_TSV = PROC / '2la_recombination_summary.tsv'
MAP_TSV = PROC / '2la_recombination_map.tsv'
RECOMB_SUMMARY = RESULTS / 'recombination_summary.tsv'
REPORT = RESULTS / 'recombination_report.md'
CAPTION = RESULTS / 'recombination_figure_caption.txt'
MANIFEST = RESULTS / 'recombination_manifest.json'
FIG = FIGURES / 'anopheles_2la_recombination_genealogy.pdf'
FIGPNG = FIGURES / 'anopheles_2la_recombination_genealogy.png'
FIG_ONLY = FIGURES / 'anopheles_2la_recombination_only.pdf'
FIG_ONLY_PNG = FIGURES / 'anopheles_2la_recombination_only.png'


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
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k, '')) for k in fields})


def sha256(path: Path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1<<20), b''):
            h.update(b)
    return h.hexdigest()


def parse_rate(value: str):
    s=value.strip().replace('~','').replace('<','')
    return float(s)


def interval_rate(marker_left, marker_right, physical_start_bp, physical_end_bp, genetic_distance_cm, karyotype_context, source):
    span_mb=(physical_end_bp-physical_start_bp)/1_000_000
    if span_mb <= 0:
        raise ValueError('physical marker interval must have positive span')
    return {
        'marker_left': marker_left,
        'marker_right': marker_right,
        'physical_start_bp': physical_start_bp,
        'physical_end_bp': physical_end_bp,
        'physical_midpoint_bp': (physical_start_bp+physical_end_bp)/2,
        'physical_span_mb': span_mb,
        'genetic_distance_cm': genetic_distance_cm,
        'recombination_cm_per_mb': genetic_distance_cm/span_mb,
        'karyotype_context': karyotype_context,
        'source': source,
    }


def classify_region(start_bp, end_bp):
    if start_bp >= INV_START and end_bp <= INV_END:
        return 'inside_2La'
    if end_bp < INV_START:
        return 'outside_left_flank'
    if start_bp > INV_END:
        return 'outside_right_flank'
    return 'boundary_overlap'


def make_summary_track():
    rows = [
        {
            'region': 'outside_left_flank',
            'start_bp': 1,
            'end_bp': INV_START-1,
            'karyotype_context': '2La/+a heterokaryotype',
            'reported_rate': '<1.0',
            'rate_relation': '<',
            'unit': 'cM/Mb',
            'source': SOURCE_CITATION,
            'notes': 'Published regional estimate for flanking regions in heterokaryotypes; not marker-level reconstruction.',
        },
        {
            'region': 'inside_2La',
            'start_bp': INV_START,
            'end_bp': INV_END,
            'karyotype_context': '2La/+a heterokaryotype',
            'reported_rate': '<0.5',
            'rate_relation': '<',
            'unit': 'cM/Mb',
            'source': SOURCE_CITATION,
            'notes': 'Published regional estimate inside rearranged region in heterokaryotypes; aligned to frozen Stage-2R 2La interval.',
        },
        {
            'region': 'outside_right_flank',
            'start_bp': INV_END+1,
            'end_bp': 49_364_325,
            'karyotype_context': '2La/+a heterokaryotype',
            'reported_rate': '<1.0',
            'rate_relation': '<',
            'unit': 'cM/Mb',
            'source': SOURCE_CITATION,
            'notes': 'Published regional estimate for flanking regions in heterokaryotypes; not marker-level reconstruction.',
        },
        {
            'region': 'chromosome_2L_background',
            'start_bp': 1,
            'end_bp': 49_364_325,
            'karyotype_context': '2L+a/2L+a homokaryotype',
            'reported_rate': '~2.0',
            'rate_relation': '~',
            'unit': 'cM/Mb',
            'source': SOURCE_CITATION,
            'notes': 'Published approximate uniform homokaryotype recombination rate on 2L; shown as contextual background only.',
        },
    ]
    write_tsv(SUMMARY_TSV, rows, ['region','start_bp','end_bp','karyotype_context','reported_rate','rate_relation','unit','source','notes'])
    return rows


def write_audits():
    RAW.mkdir(parents=True, exist_ok=True)
    SOURCE_AUDIT.write_text(f'''# Recombination source audit — Anopheles 2La

Primary source: {SOURCE_CITATION}

DOI: `{SOURCE_DOI}`

Access date: 2026-10-04

## Search order and result

1. Existing repository files: no Stump et al. marker-level recombination table was present before this extension.
2. Article DOI/Wiley landing metadata: article metadata and abstract-level recombination summaries were available; no machine-readable marker-level supplement was identified during this audit.
3. Public article records: PubMed, Sapienza IRIS, Northern Arizona University Experts, OpenAlex, and journal index pages repeat the article citation and regional summary but do not provide a numerical marker table.
4. Archival data search: no public archival marker table containing both physical and genetic positions was found during this task.

Because marker-level physical and genetic positions were not reproducibly available, this extension uses the conservative fallback: a categorical published regional-summary recombination track rather than a digitized or interpolated curve.

## Source records

| source | URL or DOI | file name | download date | checksum | data type | marker-level physical and genetic positions available |
|---|---|---|---|---|---|---|
| Stump et al. 2007 DOI | https://doi.org/{SOURCE_DOI} | none downloaded | 2026-10-04 | not applicable | article metadata/abstract summary | no |
| PubMed | https://pubmed.ncbi.nlm.nih.gov/18092999/ | none downloaded | 2026-10-04 | not applicable | article metadata/abstract summary | no |
| Sapienza IRIS | https://iris.uniroma1.it/handle/11573/363495 | none downloaded | 2026-10-04 | not applicable | institutional article record | no |
| Northern Arizona University Experts | https://experts.nau.edu/en/publications/genetic-exchange-in-2la-inversion-heterokaryotypes-of-anopheles-g/ | none downloaded | 2026-10-04 | not applicable | institutional article record | no |

No unrelated sequencing data were downloaded.
''')
    COORD_AUDIT.write_text(f'''# Recombination coordinate audit — Anopheles 2La

Current MSRC Anopheles coordinate system: `{ASSEMBLY}`.

Frozen analysis interval: `{CHROM}:{INV_START}-{INV_END}` (1-based inclusive), equivalent to 20.524058-42.165532 Mb on 2L.

Stump et al. (2007) reports experimental recombination estimates from backcross progeny using markers on 2L, but this audit did not recover a machine-readable table with marker-level physical positions and genetic-map distances. Therefore no coordinate conversion, lift-over, or marker remapping was performed.

The fallback recombination track uses the frozen MSRC 2La interval for visual alignment and labels the values as published regional estimates, not as a de novo AgamP4 marker-level recombination map. This avoids approximate coordinate rescaling from older marker/assembly systems.
''')


def write_recombination_summaries(rows):
    out=[]
    for r in rows:
        out.append({
            'mode': 'published regional summary',
            'region': r['region'],
            'karyotype_context': r['karyotype_context'],
            'n_intervals': 1,
            'mean_cm_per_mb': parse_rate(r['reported_rate']),
            'median_cm_per_mb': parse_rate(r['reported_rate']),
            'rate_relation': r['rate_relation'],
            'unit': r['unit'],
            'notes': r['notes'],
        })
    write_tsv(RECOMB_SUMMARY, out, ['mode','region','karyotype_context','n_intervals','mean_cm_per_mb','median_cm_per_mb','rate_relation','unit','notes'])
    # Marker-level file is intentionally not created when source data are unavailable.
    if MAP_TSV.exists():
        MAP_TSV.unlink()
    return out


def load_stage2_signal(path, key):
    rows=read_tsv(path)
    by={r['window_id']:r for r in rows}
    grid=read_tsv(GRID)
    out=[]
    for w in grid:
        r=by.get(w['window_id'])
        if not r: continue
        val=r.get(key, 'nan')
        try: y=float(val)
        except Exception: y=math.nan
        out.append({'window_id':w['window_id'], 'start':int(w['start']), 'end':int(w['end']), 'mid_mb':(int(w['start'])+int(w['end']))/2/1e6, 'region_class':w['region_class'], 'value':y})
    return out


def plot_recombination_panel(ax, rows, show_xlabel=False):
    ax.axvspan(INV_START/1e6, INV_END/1e6, color='#756bb1', alpha=0.16, lw=0)
    colors={'2La/+a heterokaryotype':'#d95f02', '2L+a/2L+a homokaryotype':'#1b9e77'}
    for r in rows:
        x0=int(r['start_bp'])/1e6; x1=int(r['end_bp'])/1e6
        y=parse_rate(r['reported_rate'])
        ctx=r['karyotype_context']
        if r['region']=='chromosome_2L_background':
            ax.hlines(y, x0, x1, color=colors[ctx], lw=1.8, linestyles='--', label='homokaryotype background ~2.0 cM/Mb')
        else:
            ax.hlines(y, x0, x1, color=colors[ctx], lw=5, label='heterokaryotype regional upper bound' if r['region']=='outside_left_flank' else None)
            ax.text((x0+x1)/2, y+0.07, r['reported_rate'], ha='center', va='bottom', fontsize=8, color=colors[ctx])
    ax.set_ylabel('Published recomb.\n(cM/Mb)')
    ax.set_ylim(0, 2.4)
    ax.set_xlim(0, 49.364325)
    ax.text((INV_START+INV_END)/2/1e6, 2.28, '2La', ha='center', va='top', fontsize=9, color='#4b3f8f')
    ax.text(0.01,0.92,'Stump et al. 2007 crossing estimates\nregional summary, not marker-level map', transform=ax.transAxes, va='top', fontsize=8.5, bbox=dict(facecolor='white', edgecolor='0.85', alpha=0.9))
    handles, labels=ax.get_legend_handles_labels()
    seen={}
    for h,l in zip(handles,labels):
        if l and l not in seen: seen[l]=h
    ax.legend(seen.values(), seen.keys(), loc='upper right', fontsize=8, frameon=False)
    if show_xlabel:
        ax.set_xlabel('Genomic position on 2L (Mb)')


def plot_signal(ax, signal, ylabel):
    ax.axvspan(INV_START/1e6, INV_END/1e6, color='#756bb1', alpha=0.16, lw=0)
    ax.axhline(0, color='0.45', lw=0.8, ls='--')
    xs=[r['mid_mb'] for r in signal if not math.isnan(r['value'])]
    ys=[r['value'] for r in signal if not math.isnan(r['value'])]
    colors=['#d95f02' if r['region_class']=='inside' else '#757575' if r['region_class']=='boundary' else '#1b9e77' for r in signal if not math.isnan(r['value'])]
    markers=['s' if r['region_class']=='boundary' else 'o' for r in signal if not math.isnan(r['value'])]
    for mk in ['o','s']:
        sel=[i for i,m in enumerate(markers) if m==mk]
        ax.scatter([xs[i] for i in sel], [ys[i] for i in sel], s=7 if mk=='o' else 16, c=[colors[i] for i in sel], marker=mk, linewidths=0, alpha=0.9)
    ax.set_ylabel(ylabel)
    ax.set_xlim(0,49.364325)


def make_figures(summary_rows):
    c=load_stage2_signal(CROSS, 'C')
    m=load_stage2_signal(MARGIN, 'mean_M')
    d=load_stage2_signal(TOPO, 'D')
    fig, axes=plt.subplots(4,1,figsize=(10,8.8),sharex=True,gridspec_kw={'height_ratios':[0.9,1,1,1]})
    plot_recombination_panel(axes[0], summary_rows)
    plot_signal(axes[1], c, 'C(w)')
    plot_signal(axes[2], m, 'M(w)')
    plot_signal(axes[3], d, r'D(w) = q$_A$ - q$_S$')
    axes[3].set_xlabel('Genomic position on 2L (Mb)')
    axes[0].set_title('Independent recombination-suppression evidence aligned with Anopheles 2La genealogy statistics', fontsize=12)
    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG)
    fig.savefig(FIGPNG, dpi=300)
    plt.close(fig)

    fig, ax=plt.subplots(figsize=(10,2.6))
    plot_recombination_panel(ax, summary_rows, show_xlabel=True)
    ax.set_title('Published recombination suppression across Anopheles 2La', fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG_ONLY)
    fig.savefig(FIG_ONLY_PNG, dpi=300)
    plt.close(fig)


def write_report_and_caption():
    REPORT.write_text(f'''# Recombination-suppression validation for Anopheles 2La

This supplementary analysis adds an independent recombination-suppression track to the frozen Anopheles 2La empirical results. It does not modify Stage 2R, Stage 3R, or Stage 4R biological outputs.

## Source

Primary source: {SOURCE_CITATION}

DOI: `{SOURCE_DOI}`

The source is an independent experimental crossing study. It estimated recombination from backcross progeny of `2La/+a` heterokaryotypes and `2L+a/2L+a` homokaryotype controls.

## Data availability outcome

A reproducible marker-level table containing both physical positions and genetic distances was not recovered during this audit. Therefore this extension uses a coarse published regional-summary track rather than a reconstructed marker-level recombination map. No figure digitization or interpolation was performed.

## Coordinate system

The plotted recombination regions are aligned to the frozen MSRC interval `{CHROM}:{INV_START}-{INV_END}` in the same AgamP4-style coordinate system used by Stage 2R-4R. No coordinate conversion was attempted because no marker-level source table was recovered.

## Recombination pattern

Published regional estimates report approximately `~2.0 cM/Mb` in `2L+a/2L+a` homokaryotypes, `<0.5 cM/Mb` inside 2La in `2La/+a` heterokaryotypes, and `<1.0 cM/Mb` in heterokaryotype flanking regions.

## Relationship to Stage 2R genealogy statistics

The integrated figure aligns this independent recombination-suppression evidence with the frozen Stage 2R `C(w)`, `M(w)`, and `D(w)` tracks across 2L. The intended interpretation is spatial concordance and independent biological support for the recombination-suppression component of the MSRC interpretation. It is not a causal test showing that the Stump et al. recombination estimates caused the Fontaine/MalariaGEN genealogy signal.

## Limitations

- The recombination panel is a published regional-summary track, not a marker-level recombination map.
- The recombination crossing experiment and the Fontaine/MalariaGEN genealogy analysis use different samples and study designs.
- No LD-based recombination estimator was run on the 72 Fontaine samples.
''')
    CAPTION.write_text(f'''Independent recombination-suppression evidence aligned with Anopheles 2La genealogy statistics. Panel A shows published regional recombination estimates from the crossing experiment of Stump et al. (2007), independent of the Fontaine/MalariaGEN genealogy analysis. The panel is a regional summary track, not a de novo marker-level recombination map: homokaryotype 2L+a/2L+a background recombination is approximately 2.0 cM/Mb, whereas 2La/+a heterokaryotype recombination is reported as <0.5 cM/Mb inside the inversion and <1.0 cM/Mb in flanking regions. Panels B-D show the frozen Stage 2R MalariaGEN statistics C(w), M(w), and D(w) across 50-kb windows on 2L. Purple shading marks the frozen MSRC 2La interval, 2L:{INV_START}-{INV_END}. The datasets are independent, so the figure supports spatial concordance with recombination suppression rather than a direct causal estimate.
''')


def update_readme():
    block=f'''

## Recombination-suppression validation

A supplementary recombination-suppression validation has been added using the independent crossing experiment of Stump et al. 2007 (`{SOURCE_DOI}`). No marker-level table with physical and genetic positions was reproducibly recovered, so the analysis uses a conservative published regional-summary track rather than an interpolated recombination map.

The integrated figure aligns the published recombination pattern with the frozen Stage 2R `C(w)`, `M(w)`, and `D(w)` tracks. It supports the biological assumption that recombination between alternative 2La arrangements is strongly reduced inside the inversion, while leaving Stage 2R-4R results unchanged.

Key outputs:

- `empirical/anopheles_2la/figures/anopheles_2la_recombination_genealogy.pdf`
- `empirical/anopheles_2la/results/recombination_report.md`
- `data/anopheles_2la/processed/2la_recombination_summary.tsv`
'''
    old=README.read_text() if README.exists() else ''
    marker='## Recombination-suppression validation'
    if marker in old:
        old=old.split(marker)[0].rstrip()+'\n'
    README.write_text(old.rstrip()+block)


def write_manifest():
    outputs=[SOURCE_AUDIT,COORD_AUDIT,SUMMARY_TSV,RECOMB_SUMMARY,REPORT,CAPTION,MANIFEST,FIG,FIGPNG,FIG_ONLY,FIG_ONLY_PNG,README]
    scripts=[Path(__file__)]
    stage_manifests={}
    for p in [STAGE2_MANIFEST, STAGE3_MANIFEST, STAGE4_MANIFEST]:
        if p.exists(): stage_manifests[str(p.relative_to(REPO))]=sha256(p)
    data_files={str(p.relative_to(REPO)):sha256(p) for p in [SUMMARY_TSV, RECOMB_SUMMARY] if p.exists()}
    MANIFEST.write_text(json.dumps({
        'analysis': 'Anopheles 2La recombination-suppression validation',
        'analysis_date': '2026-10-04',
        'source_citation': SOURCE_CITATION,
        'doi': SOURCE_DOI,
        'source_urls': SOURCE_URLS,
        'source_files': [],
        'source_file_hashes': {},
        'data_mode': 'published regional summary fallback; no marker-level table recovered',
        'coordinate_assembly': 'published marker-level coordinates unavailable; regional values aligned to frozen MSRC interval only',
        'current_analysis_coordinate_assembly': ASSEMBLY,
        'conversion_method': 'none; no marker-level coordinate conversion performed',
        '2La_boundaries': {'chrom':CHROM, 'start':INV_START, 'end':INV_END, 'coordinate_type':'1-based inclusive'},
        'stage2_stage3_stage4_manifest_hashes_after_extension': stage_manifests,
        'script_hashes': {str(p.relative_to(REPO)):sha256(p) for p in scripts},
        'output_hashes': {str(p.relative_to(REPO)):sha256(p) for p in outputs if p.exists()},
        'data_file_hashes': data_files,
    }, indent=2, sort_keys=True)+'\n')


def run():
    assert INV_START == 20_524_058 and INV_END == 42_165_532
    before={p:sha256(p) for p in [STAGE2_MANIFEST,STAGE3_MANIFEST,STAGE4_MANIFEST] if p.exists()}
    write_audits()
    rows=make_summary_track()
    write_recombination_summaries(rows)
    make_figures(rows)
    write_report_and_caption()
    update_readme()
    write_manifest()
    after={p:sha256(p) for p in before}
    assert before == after, 'Stage 2R/3R/4R manifests changed unexpectedly'
    print(json.dumps({'mode':'published regional summary fallback','marker_level_data_found':False,'figure':str(FIG),'summary':str(SUMMARY_TSV)}, sort_keys=True))


class Tests(unittest.TestCase):
    def test_boundaries(self):
        self.assertEqual(INV_START,20_524_058)
        self.assertEqual(INV_END,42_165_532)
    def test_classification(self):
        self.assertEqual(classify_region(INV_START,INV_END),'inside_2La')
        self.assertEqual(classify_region(1,INV_START-1),'outside_left_flank')
        self.assertEqual(classify_region(INV_END+1,49_364_325),'outside_right_flank')
        self.assertEqual(classify_region(INV_START-10,INV_START+10),'boundary_overlap')
    def test_interval_rate(self):
        r=interval_rate('a','b',1_000_000,2_000_000,2.0,'ctx','src')
        self.assertAlmostEqual(r['recombination_cm_per_mb'],2.0)
    def test_interval_rate_zero_span(self):
        with self.assertRaises(ValueError):
            interval_rate('a','b',1_000_000,1_000_000,1.0,'ctx','src')
    def test_summary_mode(self):
        rows=make_summary_track()
        self.assertEqual(len(rows),4)
        inside=[r for r in rows if r['region']=='inside_2La'][0]
        self.assertEqual(inside['reported_rate'],'<0.5')
    def test_sorted_stage2_inputs_if_present(self):
        if GRID.exists():
            mids=[]
            for r in read_tsv(GRID):
                mids.append((int(r['start'])+int(r['end']))/2)
            self.assertEqual(mids, sorted(mids))
    def test_stage_manifests_present_unchanged_during_tests(self):
        for p in [STAGE2_MANIFEST,STAGE3_MANIFEST,STAGE4_MANIFEST]:
            if p.exists():
                self.assertTrue(sha256(p))


def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument('--run-tests', action='store_true')
    args=ap.parse_args(argv)
    if args.run_tests:
        result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
        return 0 if result.wasSuccessful() else 1
    run()
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
