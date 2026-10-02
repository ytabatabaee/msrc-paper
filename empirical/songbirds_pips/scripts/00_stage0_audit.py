#!/usr/bin/env python3
"""Build Stage-0 provenance and feasibility tables for songbird PIPs.

This script intentionally uses only publication metadata, supplementary tables,
and author-reported data-availability statements. It does not read local-tree
outputs or infer any genealogy.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "songbirds_pips"
EMP = ROOT / "empirical" / "songbirds_pips"
RAW = DATA / "raw" / "oup_supplement"
XLSX = RAW / "DatasetS1_S4_S5.xlsx"
SUPP_ZIP = ROOT / "data" / "songbirds_pips_supplementary_data.zip"
TODAY = "2026-10-02"

NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


SPECIES_TABLE = [
    ("Piciformes", "Picidae", "Picoides arcticus", "Black-backed Woodpecker", 29, 0),
    ("Piciformes", "Picidae", "Dryobates villosus", "Hairy Woodpecker", 40, 0),
    ("Piciformes", "Picidae", "Sphyrapicus varius", "Yellow-bellied Sapsucker", 60, 7),
    ("Passeriformes", "Tyrannidae", "Empidonax alnorum", "Alder Flycatcher", 54, 3),
    ("Passeriformes", "Tyrannidae", "Empidonax flaviventris", "Yellow-bellied Flycatcher", 56, 5),
    ("Passeriformes", "Tyrannidae", "Empidonax minimus", "Least Flycatcher", 43, 11),
    ("Passeriformes", "Corvidae", "Perisoreus canadensis", "Canada Jay", 25, 0),
    ("Passeriformes", "Vireonidae", "Vireo olivaceus", "Red-eyed Vireo", 62, 4),
    ("Passeriformes", "Vireonidae", "Vireo philadelphicus", "Philadelphia Vireo", 15, 0),
    ("Passeriformes", "Vireonidae", "Vireo solitarius", "Blue-headed Vireo", 54, 2),
    ("Passeriformes", "Paridae", "Poecile atricapillus", "Black-capped Chickadee", 57, 2),
    ("Passeriformes", "Paridae", "Poecile hudsonicus", "Boreal Chickadee", 35, 0),
    ("Passeriformes", "Regulidae", "Corthylio calendula", "Ruby-crowned Kinglet", 47, 5),
    ("Passeriformes", "Regulidae", "Regulus satrapa", "Golden-crowned Kinglet", 59, 0),
    ("Passeriformes", "Certhiidae", "Certhia americana", "Brown Creeper", 46, 0),
    ("Passeriformes", "Troglodytidae", "Troglodytes hiemalis", "Winter Wren", 23, 2),
    ("Passeriformes", "Turdidae", "Catharus fuscescens", "Veery", 44, 18),
    ("Passeriformes", "Turdidae", "Catharus guttatus", "Hermit Thrush", 51, 8),
    ("Passeriformes", "Turdidae", "Catharus ustulatus", "Swainson's Thrush", 64, 8),
    ("Passeriformes", "Passerellidae", "Junco hyemalis", "Dark-eyed Junco", 57, 3),
    ("Passeriformes", "Passerellidae", "Melospiza lincolnii", "Lincoln's Sparrow", 51, 6),
    ("Passeriformes", "Passerellidae", "Zonotrichia albicollis", "White-throated Sparrow", 66, 6),
    ("Passeriformes", "Parulidae", "Cardellina canadensis", "Canada Warbler", 29, 7),
    ("Passeriformes", "Parulidae", "Geothlypis philadelphia", "Mourning Warbler", 50, 7),
    ("Passeriformes", "Parulidae", "Oporornis agilis", "Connecticut Warbler", 28, 4),
    ("Passeriformes", "Parulidae", "Leiothlypis peregrina", "Tennessee Warbler", 45, 7),
    ("Passeriformes", "Parulidae", "Leiothlypis ruficapilla", "Nashville Warbler", 55, 11),
    ("Passeriformes", "Parulidae", "Setophaga castanea", "Bay-breasted Warbler", 44, 10),
    ("Passeriformes", "Parulidae", "Setophaga coronata", "Yellow-rumped Warbler", 68, 8),
    ("Passeriformes", "Parulidae", "Setophaga fusca", "Blackburnian Warbler", 51, 6),
    ("Passeriformes", "Parulidae", "Setophaga magnolia", "Magnolia Warbler", 54, 5),
    ("Passeriformes", "Parulidae", "Setophaga palmarum", "Palm Warbler", 49, 5),
    ("Passeriformes", "Parulidae", "Setophaga pensylvanica", "Chestnut-sided Warbler", 46, 6),
    ("Passeriformes", "Parulidae", "Setophaga tigrina", "Cape May Warbler", 42, 3),
    ("Passeriformes", "Parulidae", "Setophaga virens", "Black-throated Green Warbler", 61, 5),
]


AUTHOR_TRANS_SPECIFIC = {
    "CM020367.1",
    "CM020368.1",
    "CM020371.1",
}


SELECTED_IDS = [
    ("Empidonax flaviventris", "CM020533.1", "A_strong_within_species"),
    ("Empidonax minimus", "CM020533.1", "B_independent_within_species"),
    ("Geothlypis philadelphia", "CM019906.1", "independent_warbler_within_species"),
    ("Zonotrichia albicollis", "CM042600.1", "D_small_LD_supported_control"),
]

SELECTED_SHARED = [
    ("shared_CM027530.1", "CM027530.1", "C_cross_species_warbler_shared_region"),
]


def ensure_dirs() -> None:
    for p in [
        DATA / "metadata",
        DATA / "processed",
        DATA / "raw" / "oup_supplement",
        EMP / "results",
        EMP / "config",
        EMP / "scripts",
        EMP / "figures",
    ]:
        p.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def xlsx_sheets(path: Path) -> dict[str, list[dict[str, str]]]:
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in ss.findall("a:si", NS):
                strings.append("".join(t.text or "" for t in si.iter(f"{{{NS['a']}}}t")))

        def cell_value(c: ET.Element) -> str:
            t = c.attrib.get("t")
            v = c.find("a:v", NS)
            if t == "inlineStr":
                is_ = c.find("a:is", NS)
                return "" if is_ is None else "".join(x.text or "" for x in is_.iter(f"{{{NS['a']}}}t"))
            if v is None:
                return ""
            val = v.text or ""
            return strings[int(val)] if t == "s" else val

        out = {}
        for s in wb.find("a:sheets", NS):
            name = s.attrib["name"]
            rid = s.attrib[f"{{{NS['r']}}}id"]
            target = rid_to_target[rid]
            if not target.startswith("worksheets/"):
                continue
            ws = ET.fromstring(z.read("xl/" + target))
            rows = []
            for row in ws.findall(".//a:sheetData/a:row", NS):
                vals = []
                last_col = 0
                for c in row.findall("a:c", NS):
                    ref = c.attrib.get("r", "")
                    letters = re.sub(r"[^A-Z]", "", ref)
                    col = 0
                    for ch in letters:
                        col = col * 26 + (ord(ch) - 64)
                    while last_col + 1 < col:
                        vals.append("")
                        last_col += 1
                    vals.append(cell_value(c))
                    last_col = col
                rows.append(vals)
            if not rows:
                out[name] = []
                continue
            header = [str(x).strip() for x in rows[0]]
            data = []
            for vals in rows[1:]:
                if not any(vals):
                    continue
                data.append({header[i]: vals[i] if i < len(vals) else "" for i in range(len(header))})
            out[name] = data
        return out


def num(x: str) -> float | None:
    if x in ("", "NA", None):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def intish(x: str) -> str:
    n = num(x)
    if n is None:
        return x or "NA"
    if abs(n - round(n)) < 1e-6:
        return str(int(round(n)))
    return str(n)


def pip_key(species: str, pip_name: str) -> str:
    return f"{species.replace(' ', '_')}|{pip_name}"


def main() -> None:
    ensure_dirs()
    sheets = xlsx_sheets(XLSX)
    d1 = sheets["DatasetS1"]
    d4 = sheets["DatasetS4"]
    d5 = sheets["DatasetS5"]

    counts = defaultdict(Counter)
    sample_rows_by_species = defaultdict(set)
    sample_rows = {}
    for r in d4:
        key = pip_key(r["Species"], r["PIP.name"])
        counts[key][r["PIP.genotype"]] += 1
        sid = r["Specimen.number"]
        sample_rows_by_species[r["Species"]].add(sid)
        sample_rows[sid] = r

    pips_by_species = defaultdict(list)
    for r in d1:
        pips_by_species[r["Species"]].append(r)

    source_rows = [
        {
            "source": "primary_article",
            "description": "Pegan and Winger 2025, Large Inversion Polymorphisms are Widespread in North American Songbirds, Genome Biology and Evolution 17(11):evaf205.",
            "url": "https://doi.org/10.1093/gbe/evaf205",
            "version_or_commit": "published 2025-11-03; corrected/typeset 2025-11-19",
            "local_path": "",
            "sha256": "",
            "notes": f"Accessed {TODAY}; article page records DOI 10.1093/gbe/evaf205.",
        },
        {
            "source": "pmc_full_text",
            "description": "PubMed Central full text for the primary article.",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12628760/",
            "version_or_commit": "PMCID PMC12628760",
            "local_path": "",
            "sha256": "",
            "notes": f"Accessed {TODAY}; PubMed record PMID 41178791.",
        },
        {
            "source": "pubmed",
            "description": "PubMed bibliographic record.",
            "url": "https://pubmed.ncbi.nlm.nih.gov/41178791/",
            "version_or_commit": "PMID 41178791",
            "local_path": "",
            "sha256": "",
            "notes": f"Accessed {TODAY}.",
        },
        {
            "source": "oup_supplement_zip",
            "description": "OUP supplementary-data ZIP containing SI PDF, DatasetS1_S4_S5.xlsx, and analysis-code text files.",
            "url": "https://academic.oup.com/gbe/article/17/11/evaf205/8313116#supplementary-data",
            "version_or_commit": "OUP supplement downloaded 2026-10-02; ZIP internal dates 2025-07-24 to 2025-11-19",
            "local_path": str(SUPP_ZIP.relative_to(ROOT)),
            "sha256": sha256(SUPP_ZIP),
            "notes": "Downloaded from the article's signed OUP CDN supplementary-data link; the stable landing page is recorded here.",
        },
        {
            "source": "supplement_workbook",
            "description": "Supplementary workbook with Dataset S1, S4, and S5.",
            "url": "https://academic.oup.com/gbe/article/17/11/evaf205/8313116#supplementary-data",
            "version_or_commit": "DatasetS1_S4_S5.xlsx from OUP supplement",
            "local_path": str(XLSX.relative_to(ROOT)),
            "sha256": sha256(XLSX),
            "notes": "Primary machine-readable source for PIP inventory, individual genotype metadata, and shared-PIP tree representatives.",
        },
        {
            "source": "supplement_pdf",
            "description": "Supplementary Information PDF.",
            "url": "https://academic.oup.com/gbe/article/17/11/evaf205/8313116#supplementary-data",
            "version_or_commit": "SI_GBE_6Oct25.pdf",
            "local_path": str((RAW / "SI_GBE_6Oct25.pdf").relative_to(ROOT)),
            "sha256": sha256(RAW / "SI_GBE_6Oct25.pdf"),
            "notes": "Contains supplementary methods/tables/figures; no tree inference was inspected for candidate selection.",
        },
        {
            "source": "authors_supplement_code",
            "description": "Analysis scripts packaged in the OUP supplement.",
            "url": "https://academic.oup.com/gbe/article/17/11/evaf205/8313116#supplementary-data",
            "version_or_commit": "packaged code, not a moving branch",
            "local_path": str((RAW / "code" / "code").relative_to(ROOT)),
            "sha256": "",
            "notes": "README states scripts are documentation of commands, not ready-to-run pipelines.",
        },
        {
            "source": "sra_bioproject_PRJNA1043688",
            "description": "NCBI SRA BioProject listed by the article for sequence reads.",
            "url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1043688",
            "version_or_commit": "NCBI BioProject PRJNA1043688",
            "local_path": "",
            "sha256": "",
            "notes": f"Not downloaded in Stage 0; article Data Availability lists this BioProject. Accessed {TODAY}.",
        },
        {
            "source": "sra_bioproject_PRJNA1130443",
            "description": "NCBI SRA BioProject listed by the article for sequence reads.",
            "url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1130443",
            "version_or_commit": "NCBI BioProject PRJNA1130443",
            "local_path": "",
            "sha256": "",
            "notes": f"Not downloaded in Stage 0; NCBI search reports raw reads, 1,897 SRA experiments, 1,867 BioSamples, and 3.37 TB. Accessed {TODAY}.",
        },
        {
            "source": "preprocessing_figshare",
            "description": "Pegan et al. 2025 preprocessing scripts cited by the supplement README.",
            "url": "https://doi.org/10.6084/m9.figshare.27284553",
            "version_or_commit": "DOI-versioned Figshare record",
            "local_path": "",
            "sha256": "",
            "notes": "The songbird PIP supplement README points to this Figshare DOI for raw-read preprocessing scripts.",
        },
        {
            "source": "local_pcangsd",
            "description": "External local PCAngsd code referenced by authors.",
            "url": "https://github.com/alxsimon/local_pcangsd",
            "version_or_commit": "not frozen locally",
            "local_path": "",
            "sha256": "",
            "notes": "Authors adapted this tool; Stage 0 depends on the packaged author outputs, not the moving GitHub branch.",
        },
    ]
    write_tsv(DATA / "metadata" / "source_manifest.tsv", source_rows, ["source", "description", "url", "version_or_commit", "local_path", "sha256", "notes"])

    pip_inventory = []
    for r in d1:
        key = pip_key(r["Species"], r["PIP.name"])
        c = counts[key]
        pca_n = int(float(r["PCA.cluster.number"])) if r["PCA.cluster.number"] not in ("", "NA") else ""
        ld = r["LD.supported"] == "1"
        shared = "yes" if "multispecies" in r["PIP.unique.status"] else "no"
        author_conf = "PCA_plus_LD" if ld else "PCA_only_LD_ambiguous"
        pip_inventory.append({
            "pip_id": key,
            "species": r["Species"],
            "scientific_name": r["Species"],
            "reference_genome_assembly": r["Reference.genome.assembly"],
            "reference_genome_species": r["Reference.genome.species"],
            "chromosome": r["Chromosome.accession"],
            "start": intish(r["PIP.start"]),
            "end": intish(r["PIP.end"]),
            "length_bp": str(int(round((num(r["PIP.length"]) or 0) * 1_000_000))),
            "length_mb": r["PIP.length"],
            "pca_support": "yes",
            "ld_support": "yes" if ld else "ambiguous",
            "other_support": f"PCA_clusters={pca_n}; genotype_counts={dict(c)}; IBD_slope={r['IBD.slope'] or 'NA'}; min_HOR_dist={r['min_HOR_dist'] or 'NA'}",
            "author_confidence": author_conf,
            "shared_across_species": shared,
            "notes": f"authors_PIP_unique_status={r['PIP.unique.status']}; comp_break_right={intish(r['comp_break_right'])}; comp_break_left={intish(r['comp_break_left'])}",
        })
    write_tsv(
        DATA / "processed" / "stage0_pip_inventory.tsv",
        pip_inventory,
        ["pip_id", "species", "scientific_name", "reference_genome_assembly", "reference_genome_species", "chromosome", "start", "end", "length_bp", "length_mb", "pca_support", "ld_support", "other_support", "author_confidence", "shared_across_species", "notes"],
    )

    genotype_rows = []
    for r in d1:
        key = pip_key(r["Species"], r["PIP.name"])
        c = counts[key]
        n = sum(c.values())
        genotype_rows.append({
            "pip_id": key,
            "species": r["Species"],
            "chromosome": r["Chromosome.accession"],
            "n_total": n,
            "n_standard_hom": c.get("AA", 0),
            "n_heterozygous": c.get("AB", 0),
            "n_inverted_hom": c.get("BB", 0),
            "n_unclassified": n - c.get("AA", 0) - c.get("AB", 0) - c.get("BB", 0),
            "author_genotype_counts": ";".join(f"{k}:{v}" for k, v in sorted(c.items())),
            "mapping_notes": "AA is the more frequent homozygote; AB is the middle PCA cluster; BB is the less frequent homozygote for 3-cluster PIPs. C-H are author arbitrary labels for non-3-cluster PIPs.",
        })
    write_tsv(DATA / "processed" / "stage0_pip_genotype_counts.tsv", genotype_rows, ["pip_id", "species", "chromosome", "n_total", "n_standard_hom", "n_heterozygous", "n_inverted_hom", "n_unclassified", "author_genotype_counts", "mapping_notes"])

    sample_inventory = []
    sp_ref = {r["Species"]: (r["Reference.genome.assembly"], r["Reference.genome.species"]) for r in d1}
    for sid, r in sorted(sample_rows.items()):
        ref = sp_ref.get(r["Species"], ("", ""))
        sample_inventory.append({
            "sample_id": sid,
            "species": r["Species"],
            "population": "not_provided; latitude_longitude_available",
            "sex": "not_provided",
            "assembly_or_reference": f"{ref[0]} ({ref[1]})",
            "sequencing_source": f"NCBI SRA BioSample/accession {r['SRA.Accession']}",
            "inversion_genotype_available": "yes_for_at_least_one_PIP",
            "source": "DatasetS4",
            "voucher_museum": r["Voucher.museum"],
            "latitude": r["Latitude"],
            "longitude": r["Longitude"],
            "sra_accession": r["SRA.Accession"],
        })
    write_tsv(DATA / "processed" / "stage0_sample_inventory.tsv", sample_inventory, ["sample_id", "species", "population", "sex", "assembly_or_reference", "sequencing_source", "inversion_genotype_available", "source", "voucher_museum", "latitude", "longitude", "sra_accession"])

    species_rows = []
    for order, fam, sp, common, n, npip in SPECIES_TABLE:
        sra_count = len(sample_rows_by_species.get(sp, set()))
        species_rows.append({
            "species": sp,
            "common_name": common,
            "order": order,
            "family": fam,
            "total_individuals_article_table1": n,
            "individuals_in_datasetS4": sra_count,
            "number_of_pips_article_table1": npip,
            "inversion_genotype_assignments": "yes" if npip > 0 else "no_pips_detected; not in DatasetS4",
            "raw_reads_available": "yes_SRA_PRJNA1043688_or_PRJNA1130443",
            "bam_cram_available": "not_publicly_packaged",
            "vcf_available": "not_publicly_packaged",
            "phased_data_available": "not_publicly_packaged",
            "reference_genome_available": "yes_GenBank_accession_in_DatasetS1" if sp in sp_ref else "not_recorded_in_DatasetS1_for_no-PIP_species",
            "notes": "DatasetS4 excludes the 7 species with no detected PIPs." if npip == 0 else "",
        })
    write_tsv(DATA / "processed" / "stage0_species_inventory.tsv", species_rows, ["species", "common_name", "order", "family", "total_individuals_article_table1", "individuals_in_datasetS4", "number_of_pips_article_table1", "inversion_genotype_assignments", "raw_reads_available", "bam_cram_available", "vcf_available", "phased_data_available", "reference_genome_available", "notes"])

    shared_groups = defaultdict(list)
    for r in d5:
        shared_groups[r["Shared.PIP.tree"]].append(r)
    cross_rows = []
    for shared_id, rows in sorted(shared_groups.items()):
        species = sorted(set(r["Species"] for r in rows))
        pip_names = sorted(set(r["PIP.name"] for r in rows))
        d1_rows = [r for r in d1 if r["PIP.name"] in pip_names and r["Species"] in species]
        chroms = sorted(set(r["Chromosome.accession"] for r in d1_rows)) or [shared_id]
        refs = sorted(set(f"{r['Reference.genome.assembly']} ({r['Reference.genome.species']})" for r in d1_rows))
        cross_rows.append({
            "shared_pip_id": f"shared_{shared_id}",
            "species_involved": ";".join(species),
            "chromosome_or_region": ";".join(chroms),
            "species_level_pip_names": ";".join(pip_names),
            "evidence_for_homology": "authors_identified_overlapping_PIPs_mapped_to_same_reference_chromosome; DatasetS5 representative set available",
            "inversion_genotype_information": ";".join(f"{k}:{v}" for k, v in sorted(Counter(r["PIP.geno"] for r in rows).items())),
            "same_reference_coordinate_system": "yes_for_species_mapped_to_same_reference; verify before cross-species alignment",
            "sequence_tree_analysis_feasible": "possible_from_SRA_reads_after_recreating_BAM_SAF_distance_pipeline; no reusable Newick/data matrix released",
            "biological_interpretation_by_authors": "author-reported shared/candidate trans-specific region; do not infer introgression versus ancient polymorphism without Stage-1 tests",
            "notes": f"DatasetS5 representative_rows={len(rows)}; reference_genomes={';'.join(refs)}",
        })
    write_tsv(DATA / "processed" / "stage0_cross_species_pips.tsv", cross_rows, ["shared_pip_id", "species_involved", "chromosome_or_region", "species_level_pip_names", "evidence_for_homology", "inversion_genotype_information", "same_reference_coordinate_system", "sequence_tree_analysis_feasible", "biological_interpretation_by_authors", "notes"])

    availability_rows = [
        {"data_object": "existing local phylogenetic trees", "available": "partial_figures_only", "format": "supplementary plots/PDF, not reusable Newick", "scope": "shared-PIP NJ visualizations reported by authors", "sample_count": "representative subsets", "species_count": "varies", "size": "not machine-readable as trees", "source": "article Dataset S3 / methods", "usable_for_local_phylogeny": "no_for_rigorous_MSRC_reanalysis", "notes": "Do not use plotted topology for Stage-0 candidate selection."},
        {"data_object": "whole-genome trees", "available": "no", "format": "not found", "scope": "none", "sample_count": "", "species_count": "", "size": "", "source": "OUP supplement/code audit", "usable_for_local_phylogeny": "no", "notes": "Background species topology must be estimated or sourced independently in Stage 1."},
        {"data_object": "SNP matrices", "available": "no", "format": "not packaged", "scope": "authors generated genotype likelihood inputs", "sample_count": "1660", "species_count": "35", "size": "large", "source": "supplement README/code", "usable_for_local_phylogeny": "not_without_reprocessing", "notes": "Scripts refer to BEAGLE/genotype-likelihood inputs but do not release them."},
        {"data_object": "VCF/BCF files", "available": "no", "format": "not packaged", "scope": "none public in supplement", "sample_count": "", "species_count": "", "size": "", "source": "OUP supplement/code audit", "usable_for_local_phylogeny": "no", "notes": "No analysis-ready VCF was identified."},
        {"data_object": "genotype likelihoods", "available": "no", "format": "ANGSD/PCAngsd intermediate inputs referenced only", "scope": "per species/chromosome", "sample_count": "1660 if regenerated", "species_count": "35", "size": "large", "source": "author scripts", "usable_for_local_phylogeny": "not_without_reprocessing", "notes": "Cheaper than raw FASTQ only if authors can provide intermediate files; not public in supplement."},
        {"data_object": "phased haplotypes", "available": "no", "format": "not found", "scope": "none", "sample_count": "", "species_count": "", "size": "", "source": "OUP supplement/code audit", "usable_for_local_phylogeny": "no", "notes": "The authors used low-coverage genotype likelihoods, not phased haplotypes."},
        {"data_object": "sequence alignments", "available": "no", "format": "not packaged", "scope": "none", "sample_count": "", "species_count": "", "size": "", "source": "OUP supplement/code audit", "usable_for_local_phylogeny": "no", "notes": "Local alignments would need reconstruction from reads/BAMs."},
        {"data_object": "BAM/CRAM files", "available": "not_publicly_packaged", "format": "BAMs referenced as local inputs", "scope": "all samples after preprocessing", "sample_count": "1660", "species_count": "35", "size": "very large", "source": "supplement README and scripts", "usable_for_local_phylogeny": "yes_if_obtained_or_regenerated", "notes": "Authors' B3 pipeline requires BAM lists; no public BAM archive found in Stage 0."},
        {"data_object": "raw reads", "available": "yes", "format": "NCBI SRA FASTQ/raw reads", "scope": "sequence reads", "sample_count": "1660 used here; PRJNA1130443 has more experiments/samples", "species_count": "35", "size": "terabytes", "source": "PRJNA1043688; PRJNA1130443; DatasetS4 SRA accessions", "usable_for_local_phylogeny": "yes_but_expensive", "notes": "Reliable but costly route; avoid full reprocessing unless no intermediate data can be obtained."},
        {"data_object": "reference genomes", "available": "yes", "format": "GenBank assemblies/accessions", "scope": "reference per mapped species group", "sample_count": "", "species_count": "28 PIP species in DatasetS1", "size": "large if downloaded", "source": "DatasetS1", "usable_for_local_phylogeny": "yes", "notes": "Coordinates are reference-specific and must not be merged across assemblies without liftover/homology evidence."},
        {"data_object": "PIP coordinates", "available": "yes", "format": "xlsx/TSV", "scope": "174 PIPs", "sample_count": "", "species_count": "28 with PIPs", "size": "small", "source": "DatasetS1", "usable_for_local_phylogeny": "yes", "notes": "Approximate breakpoints mapped to reference chromosomes."},
        {"data_object": "inversion genotype assignments", "available": "yes", "format": "xlsx/TSV", "scope": "sample-by-PIP assignments", "sample_count": str(len(sample_rows)), "species_count": "28 with PIPs", "size": "small", "source": "DatasetS4", "usable_for_local_phylogeny": "yes_for_labeling_tests", "notes": "AA/AB/BB for 3-cluster PIPs; C-H arbitrary labels for non-3-cluster PIPs."},
        {"data_object": "PCA classifications", "available": "yes", "format": "genotype labels and PCA cluster counts", "scope": "174 PIPs", "sample_count": str(len(sample_rows)), "species_count": "28", "size": "small", "source": "DatasetS1/DatasetS4", "usable_for_local_phylogeny": "labels_only", "notes": "Not sufficient as phylogenetic input."},
        {"data_object": "LD evidence", "available": "yes", "format": "binary LD.supported field plus plots", "scope": "174 PIPs", "sample_count": "", "species_count": "28", "size": "small", "source": "DatasetS1/Supplement", "usable_for_local_phylogeny": "selection_metadata_only", "notes": "139 PIPs are LD-supported."},
        {"data_object": "species/population metadata", "available": "partial", "format": "sample IDs, museum, coordinates, SRA", "scope": "PIP species only in DatasetS4", "sample_count": str(len(sample_rows)), "species_count": "28", "size": "small", "source": "DatasetS4; no-PIP species metadata in Pegan et al. 2025 supplement", "usable_for_local_phylogeny": "yes_for_covariates", "notes": "Population labels and sex are not provided in DatasetS4."},
    ]
    write_tsv(EMP / "results" / "stage0_data_availability.tsv", availability_rows, ["data_object", "available", "format", "scope", "sample_count", "species_count", "size", "source", "usable_for_local_phylogeny", "notes"])

    inv_by_key = {r["pip_id"]: r for r in pip_inventory}
    geno_by_key = {r["pip_id"]: r for r in genotype_rows}
    candidates = []
    for inv in pip_inventory:
        g = geno_by_key[inv["pip_id"]]
        n_classes = len([x for x in g["author_genotype_counts"].split(";") if x])
        aa, ab, bb = int(g["n_standard_hom"]), int(g["n_heterozygous"]), int(g["n_inverted_hom"])
        good_balance = aa >= 5 and ab >= 5 and bb >= 5
        ld = inv["ld_support"] == "yes"
        pca3 = "PCA_clusters=3" in inv["other_support"]
        comp = "comp_break_right=NA" not in inv["notes"] or "comp_break_left=NA" not in inv["notes"]
        feasible = ld and pca3 and int(g["n_total"]) >= 25 and good_balance
        priority = "high" if feasible else "medium" if ld and int(g["n_total"]) >= 20 else "low"
        if not comp:
            priority = "low"
        candidates.append({
            "pip_id": inv["pip_id"],
            "species": inv["species"],
            "n_samples": g["n_total"],
            "n_genotype_classes": n_classes,
            "pca_supported": "yes",
            "ld_supported": "yes" if ld else "ambiguous",
            "sequence_data_available": "raw_SRA_reads_only_public",
            "background_data_available": "COMP_coordinates_available" if comp else "no_clear_COMP_region",
            "cross_species": inv["shared_across_species"],
            "phylogeny_feasible": "yes_but_requires_BAM_or_raw_read_reprocessing" if feasible else "limited",
            "priority_class": priority,
            "reason": f"LD={ld}; PCA3={pca3}; n={g['n_total']}; AA/AB/BB={aa}/{ab}/{bb}; comp_region={'yes' if comp else 'no'}",
        })
    for shared_id, shared_tree, label in SELECTED_SHARED:
        rows = shared_groups[shared_tree]
        candidates.append({
            "pip_id": shared_id,
            "species": ";".join(sorted(set(r["Species"] for r in rows))),
            "n_samples": len(rows),
            "n_genotype_classes": len(set(r["PIP.geno"] for r in rows)),
            "pca_supported": "yes_for_constituent_species",
            "ld_supported": "mixed_verify_constituent_species",
            "sequence_data_available": "raw_SRA_reads_only_public",
            "background_data_available": "authors_define_shared_COMP_strategy; regenerate required",
            "cross_species": "yes",
            "phylogeny_feasible": "yes_but_requires_regenerating_author_distance_pipeline",
            "priority_class": "high",
            "reason": f"{label}; DatasetS5 has representative AA/BB samples across {len(set(r['Species'] for r in rows))} species; selected before any Stage-0 tree inspection.",
        })
    write_tsv(EMP / "results" / "stage0_candidate_pips.tsv", candidates, ["pip_id", "species", "n_samples", "n_genotype_classes", "pca_supported", "ld_supported", "sequence_data_available", "background_data_available", "cross_species", "phylogeny_feasible", "priority_class", "reason"])

    selected = []
    for sp, pip, label in SELECTED_IDS:
        key = pip_key(sp, pip)
        inv, g = inv_by_key[key], geno_by_key[key]
        selected.append((key, label, inv, g))
    overview = []
    for key, label, inv, g in selected:
        overview.append({
            "PIP": key,
            "Species": inv["species"],
            "Chr": inv["chromosome"],
            "Size": f"{float(inv['length_mb']):.3g} Mb",
            "N": g["n_total"],
            "PCA": "yes",
            "LD": inv["ld_support"],
            "Cross-species": inv["shared_across_species"],
            "Phylogeny input": "SRA reads; regenerate BAM/SAF or request intermediates",
            "Priority": "high" if "control" not in label else "medium",
        })
    rows = shared_groups["CM027530.1"]
    overview.append({
        "PIP": "shared_CM027530.1",
        "Species": f"{len(set(r['Species'] for r in rows))} warbler species",
        "Chr": "CM027530.1",
        "Size": "shared overlap; derive in Stage 1",
        "N": len(rows),
        "PCA": "yes",
        "LD": "verify per species",
        "Cross-species": "yes",
        "Phylogeny input": "SRA reads; regenerate shared-PIP distance pipeline",
        "Priority": "high",
    })
    write_tsv(EMP / "results" / "stage0_dataset_summary.tsv", overview, ["PIP", "Species", "Chr", "Size", "N", "PCA", "LD", "Cross-species", "Phylogeny input", "Priority"])

    with (EMP / "results" / "stage0_dataset_summary.md").open("w") as f:
        fields = ["PIP", "Species", "Chr", "Size", "N", "PCA", "LD", "Cross-species", "Phylogeny input", "Priority"]
        f.write("| " + " | ".join(fields) + " |\n")
        f.write("|" + "|".join(["---"] * len(fields)) + "|\n")
        for r in overview:
            f.write("| " + " | ".join(str(r.get(k, "")) for k in fields) + " |\n")

    total_pips = len(d1)
    ld_count = sum(1 for r in d1 if r["LD.supported"] == "1")
    nonld_count = total_pips - ld_count
    pca3_count = sum(1 for r in d1 if r["PCA.cluster.number"] == "3")
    pca2_count = sum(1 for r in d1 if r["PCA.cluster.number"] == "2")
    pca_multi_count = total_pips - pca3_count - pca2_count
    with (EMP / "results" / "stage0_report.md").open("w") as f:
        f.write("# Songbirds PIPs Stage-0 feasibility and provenance audit\n\n")
        f.write("This Stage-0 audit freezes source provenance, PIP/sample metadata, and topology-independent candidate-selection criteria for Pegan and Winger (2025). No local trees were inferred, and no candidate was selected by inspecting Stage-0 tree outcomes.\n\n")
        f.write("## Dataset size and released data\n\n")
        f.write(f"- Primary article: Pegan and Winger, Genome Biology and Evolution 17(11):evaf205, DOI 10.1093/gbe/evaf205.\n")
        f.write(f"- Public sequence data: NCBI SRA BioProjects PRJNA1043688 and PRJNA1130443.\n")
        f.write(f"- Supplementary data downloaded locally: `data/songbirds_pips/raw/oup_supplement/`.\n")
        f.write(f"- Species examined: {len(SPECIES_TABLE)}.\n")
        f.write(f"- Individuals in article Table 1: {sum(x[4] for x in SPECIES_TABLE)}.\n")
        f.write(f"- Individuals with PIP genotype rows in Dataset S4: {len(sample_rows)} across {len(sample_rows_by_species)} PIP-containing species.\n")
        f.write(f"- Total PIPs in Dataset S1: {total_pips}; LD-supported/high-confidence PIPs: {ld_count}; PCA-only/LD-ambiguous PIPs: {nonld_count}.\n")
        f.write(f"- PCA cluster classes: {pca3_count} three-cluster PIPs, {pca2_count} two-cluster PIPs, {pca_multi_count} four-to-six-cluster PIPs.\n\n")
        f.write("## Feasibility answer\n\n")
        f.write("Local genealogies can be reconstructed, but not directly from packaged analysis-ready phylogenetic inputs. The cheapest reliable route is to request or recover the authors' intermediate BAM/ANGSD genotype-likelihood/SAF files. Without those, Stage 1 must download SRA reads for selected samples only, regenerate BAMs/genotype likelihoods for candidate regions and matched controls, and then infer local genealogies or distance trees. PCA output alone is label metadata and is not a phylogenetic input.\n\n")
        f.write("## Data availability summary\n\n")
        f.write("- Available now: PIP coordinates, reference assemblies, LD-support flags, PCA cluster counts, per-individual PIP genotype assignments, sample coordinates/SRA accessions for PIP species, and representative sample lists for shared-PIP author trees.\n")
        f.write("- Not found as public packaged analysis-ready files: VCF/BCF, phased haplotypes, SNP matrices, BAM/CRAM, alignments, reusable Newick local trees, or whole-genome trees.\n")
        f.write("- Raw reads are public but terabyte-scale; Stage 0 did not download them.\n\n")
        f.write("## Candidate-selection criteria\n\n")
        f.write("Candidates were scored using only structural/data-quality criteria: PCA support, LD support, genotype-class representation, sample count, clear coordinates, COMP/control feasibility, SRA availability, and author-reported shared status. No topology-derived statistic from local trees was used. The one cross-species candidate uses author-reported overlap/shared-PIP metadata; Stage 0 did not inspect tree topology plots.\n\n")
        f.write("Anti-circularity audit: Stage 0 read the primary article's prose summary of the authors' neighbor-joining/shared-PIP analysis while inventorying data availability and cross-species claims. No Dataset S3 tree plots, Newick trees, distance matrices, or newly inferred local trees were inspected. Candidate selection used PIP coordinates, PCA/LD/genotype metadata, sample availability, and author-declared shared-region membership only; it did not use whether a candidate's tree looked unusual or favorable for MSRC.\n\n")
        f.write("## Selected initial PIPs\n\n")
        for key, label, inv, g in selected:
            f.write(f"- `{key}` ({label}): {inv['chromosome']}:{inv['start']}-{inv['end']}, {float(inv['length_mb']):.3g} Mb, n={g['n_total']}, AA/AB/BB={g['n_standard_hom']}/{g['n_heterozygous']}/{g['n_inverted_hom']}, LD={inv['ld_support']}, cross-species={inv['shared_across_species']}.\n")
        f.write(f"- `shared_CM027530.1` (C_cross_species_warbler_shared_region): {len(set(r['Species'] for r in rows))} warbler species, DatasetS5 representative n={len(rows)}, author shared-PIP region on CM027530.1; exact shared-overlap coordinates should be reconstructed from constituent PIPs before Stage 1.\n\n")
        f.write("## Cross-species/shared PIPs\n\n")
        f.write(f"Dataset S5 contains {len(shared_groups)} shared-PIP tree groups and {len(d5)} representative rows. The article reports 24 chromosomes with PIPs in multiple species, 26 overlapping PIP cases, and 12 author-interpreted trans-specific cases. The processed cross-species inventory records all DatasetS5 shared groups without inferring introgression or ancient balanced polymorphism beyond the authors' framing.\n\n")
        f.write("## Matched control strategy\n\n")
        f.write("For within-species candidates, use the authors' COMP intervals where available, then add matched same-chromosome windows outside PIPs with comparable length/SNP density/missingness after Stage-1 genotyping. For the cross-species candidate, reconstruct the authors' shared-PIP and shared-COMP design, then add matched non-PIP regions on the same reference chromosome only after confirming all species use the same reference coordinate frame.\n\n")
        f.write("## Coordinate and assembly cautions\n\n")
        f.write("Coordinates are species/reference-specific. Dataset S1 gives the reference genome assembly and reference species for each PIP. Do not merge regions across species unless the species were mapped to the same reference genome and overlapping coordinates are explicitly reconstructed. No liftover resources were identified in Stage 0.\n\n")
        f.write("## Computational estimate\n\n")
        f.write("A selected within-species PIP would require tens of individuals and one focal interval plus matched controls; local execution may be reasonable if starting from BAM/CRAM or small regional read extractions, but full SRA-to-BAM processing is HPC-scale. The full public read archive is terabyte-scale. The cross-species candidate spans many species and should be run on HPC unless author intermediates are obtained.\n\n")
        f.write("## Blockers\n\n")
        f.write("- No public analysis-ready VCF/BAM/CRAM/phased haplotype/local-tree files were identified.\n")
        f.write("- Dataset S4 lacks sex and explicit population labels; only coordinates and museums are provided.\n")
        f.write("- Seven no-PIP species are not included in Dataset S4; their full metadata are in the related Pegan et al. 2025 materials.\n")
        f.write("- Shared-PIP exact overlap coordinates for cross-species tests need reconstruction from constituent PIPs and/or author metadata before analysis.\n\n")
        f.write("## Stage 1 recommendation\n\n")
        f.write("Do not launch genome-wide processing. First contact/check for author-provided BAM, ANGSD genotype-likelihood, SAF, or regional distance-matrix intermediates. If unavailable, build a small SRA-reprocessing pilot for the selected PIPs only, starting with one within-species candidate and one matched COMP/control region. Freeze the exact sample sets and control intervals before any local-tree inference.\n")

    readme = EMP / "README.md"
    readme.write_text("# Songbirds PIPs empirical analysis\n\nStage-0 feasibility/provenance audit for Pegan and Winger (2025) North American songbird putative inversion polymorphisms. This directory intentionally contains no tree inference yet.\n")
    (DATA / "README.md").write_text("# Songbirds PIPs data\n\nRaw supplementary data and Stage-0 processed metadata for Pegan and Winger (2025). Large SRA read datasets were not downloaded during Stage 0.\n")

    print(f"wrote {len(pip_inventory)} PIPs, {len(sample_inventory)} samples, {len(cross_rows)} shared-PIP groups")
    print(f"selected: {', '.join([x[0] for x in selected] + ['shared_CM027530.1'])}")


if __name__ == "__main__":
    main()
