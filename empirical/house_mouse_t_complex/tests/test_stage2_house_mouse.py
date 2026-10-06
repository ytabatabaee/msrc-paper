"""Lightweight Stage-2 tests that do not require the upstream archive."""

from __future__ import annotations

import sys
import importlib.util
from pathlib import Path

import pytest
from PIL import Image

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_stage2_utils import (  # noqa: E402
    Q_OTHER,
    Q_SPECIES,
    Q_T_ALT,
    Q_UNRESOLVED,
    four_taxon_topology_from_newick,
    quartet_counts_for_tree,
    spatial_thin_records,
    tips_by_species,
    tree_bipartitions,
    annotation_for_split,
    read_tsv,
)


LABELS = {
    "Mus musculus domesticus": ["dom"],
    "Mus musculus musculus": ["mus"],
    "Mus musculus castaneus": ["cast"],
    "Mus spretus": ["spret"],
}


@pytest.mark.parametrize(
    ("newick", "expected"),
    [
        ("((mus,cast),(dom,spret));", Q_SPECIES),
        ("((dom,mus),(cast,spret));", Q_T_ALT),
        ("((dom,cast),(mus,spret));", Q_OTHER),
        ("(dom,mus,cast,spret);", Q_UNRESOLVED),
    ],
)
def test_exact_quartet_classification(newick: str, expected: str) -> None:
    counts = quartet_counts_for_tree(newick, LABELS)
    assert counts[expected] == 1
    assert sum(counts.values()) == 1


def test_four_taxon_ordering_is_unrooted() -> None:
    species_newick = "((Mus_musculus_musculus,Mus_musculus_castaneus),(Mus_musculus_domesticus,Mus_spretus));"
    reordered = "((Mus_spretus,Mus_musculus_domesticus),(Mus_musculus_castaneus,Mus_musculus_musculus));"
    assert four_taxon_topology_from_newick(species_newick) == Q_SPECIES
    assert four_taxon_topology_from_newick(reordered) == Q_SPECIES


def test_population_split_comparison_ignores_child_order() -> None:
    a = tree_bipartitions("((A,B),(C,D));")
    b = tree_bipartitions("((D,C),(B,A));")
    assert a == b
    assert tree_bipartitions("((A,B),(C,D));") == tree_bipartitions("((D,C),(B,A));")


def test_treatment_filters_require_source_confidence() -> None:
    mapping = [
        {"tree_tip": "dom_std", "subspecies": "Mus musculus domesticus", "t_status": "standard_noncarrier", "mapping_confidence": "strong"},
        {"tree_tip": "dom_t", "subspecies": "Mus musculus domesticus", "t_status": "pseudo-t_haplotype", "mapping_confidence": "strong"},
        {"tree_tip": "mus_unresolved", "subspecies": "Mus musculus musculus", "t_status": "standard_noncarrier", "mapping_confidence": "tentative"},
        {"tree_tip": "spret", "subspecies": "Mus spretus", "t_status": "outgroup_not_t_haplotype", "mapping_confidence": "exact"},
    ]
    assert tips_by_species(mapping, "STANDARD_ONLY")["Mus musculus domesticus"] == ["dom_std"]
    assert tips_by_species(mapping, "T_ONLY")["Mus musculus domesticus"] == ["dom_t"]
    assert "mus_unresolved" not in tips_by_species(mapping, "STANDARD_ONLY").get("Mus musculus musculus", [])


def test_spatial_thinning_is_deterministic_and_respects_gaps() -> None:
    meta = [
        {"start_bp": "5000000"},
        {"start_bp": "5005000"},
        {"start_bp": "5200000"},
        {"start_bp": "5205000"},
    ]
    first = spatial_thin_records(meta, 10000, 0)
    second = spatial_thin_records(meta, 10000, 0)
    assert first == second
    assert first == [0, 2]


def test_main_figure_is_nonempty() -> None:
    figure = Path(__file__).resolve().parents[1] / "figures" / "house_mouse_t_complex_stage2_main.png"
    assert figure.stat().st_size > 20_000
    with Image.open(figure).convert("RGB") as image:
        pixels = list(image.getdata())
    assert image.width >= 800 and image.height >= 600
    assert len(set(pixels)) > 100
    assert sum(any(channel < 245 for channel in pixel) for pixel in pixels) > 5_000


def test_split_specific_annotation_and_normalization() -> None:
    newick = "((Mus_musculus_musculus,Mus_musculus_castaneus)'[CULength=0.4;SULength=0.2;localPP=0.9;q1=0.5;q2=0.3;q3=0.2;f1=5;f2=3;f3=2]':0.4,(Mus_musculus_domesticus,Mus_spretus));"
    target = (("Mus_musculus_castaneus", "Mus_musculus_musculus"), ("Mus_musculus_domesticus", "Mus_spretus"))
    annotation = annotation_for_split(newick, target)
    assert annotation is not None
    assert annotation["CULength"] == 0.4
    assert abs(sum(float(annotation[key]) for key in ("q1", "q2", "q3")) - 1) < 1e-8
    spec = importlib.util.spec_from_file_location("filtering", SCRIPT_DIR / "02d_filtering_robustness.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    assert module.normalize_filter_tip("X_tHaplSubsetOverallCovFiltered") == "X_tHaplSubset.fa"
    assert module.normalize_filter_tip("X_OverallCovFiltered_DelRemoved") == "X"
    assert module.normalize_filter_tip("X_typo") == "X_typo"


def test_all_eight_arrangement_patterns_and_pooling() -> None:
    patterns = {a + b + c for a in "ST" for b in "ST" for c in "ST"}
    assert patterns == {"SSS", "SST", "STS", "STT", "TSS", "TST", "TTS", "TTT"}
    assert {p for p in patterns if p.count("T") == 1} == {"SST", "STS", "TSS"}
    assert {p for p in patterns if p.count("T") == 2} == {"STT", "TST", "TTS"}


def test_arrangement_status_quartet_counts_are_exact() -> None:
    spec = importlib.util.spec_from_file_location("arrangements", SCRIPT_DIR / "02g_arrangement_status_quartets.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    mapping = [
        {"tree_tip": "dom_S", "subspecies": "Mus musculus domesticus", "t_status": "standard_noncarrier", "mapping_confidence": "strong"},
        {"tree_tip": "dom_T", "subspecies": "Mus musculus domesticus", "t_status": "pseudo-t_haplotype", "mapping_confidence": "strong"},
        {"tree_tip": "mus_S", "subspecies": "Mus musculus musculus", "t_status": "standard_noncarrier", "mapping_confidence": "strong"},
        {"tree_tip": "mus_T", "subspecies": "Mus musculus musculus", "t_status": "pseudo-t_haplotype", "mapping_confidence": "strong"},
        {"tree_tip": "cast_S", "subspecies": "Mus musculus castaneus", "t_status": "standard_noncarrier", "mapping_confidence": "strong"},
        {"tree_tip": "cast_T", "subspecies": "Mus musculus castaneus", "t_status": "pseudo-t_haplotype", "mapping_confidence": "strong"},
        {"tree_tip": "spret", "subspecies": "Mus spretus", "t_status": "outgroup_not_t_haplotype", "mapping_confidence": "exact"},
    ]
    counts = module.count_pattern("((mus_S,cast_S),(dom_S,spret));", "SSS", mapping)
    assert counts[Q_SPECIES] == 1 and sum(counts.values()) == 1


def test_full_phase_enumeration() -> None:
    assert len(range(0, 100_000, 5_000)) == 20
    assert len(range(0, 250_000, 5_000)) == 50
    assert len(range(0, 500_000, 5_000)) == 100


def test_astral_exact_validation_and_population_split_statuses() -> None:
    validation = read_tsv(Path(__file__).resolve().parents[1] / "results" / "stage2_balanced_quartet_astral_validation.tsv")
    assert len(validation) == 60
    assert all(row["concordant"] == "true" for row in validation)
    comparison = read_tsv(Path(__file__).resolve().parents[1] / "results" / "stage2_population_split_comparison.tsv")
    assert sum(row["status"] == "lost" for row in comparison) == 1
    assert sum(row["status"] == "gained" for row in comparison) == 1
    assert sum(row["status"] == "shared" for row in comparison) == 3


def test_stage2_visualization_outputs_and_numeric_invariants() -> None:
    root = Path(__file__).resolve().parents[1]
    validation = read_tsv(root / "results" / "stage2_visualization_validation.tsv")
    assert len(validation) == 7
    assert all(int(row["nonwhite_pixels"]) > 5_000 for row in validation)
    qrows = read_tsv(root / "results" / "stage2_fixed_quartet_scan.tsv")
    assert all(abs(sum(float(row[key]) for key in ("q_species", "q_t_alt", "q_other", "q_unresolved")) - 1) < 1e-8 for row in qrows)
    recomb = read_tsv(root / "results" / "stage2_recombination_state_500kb.tsv")
    assert all(abs(sum(float(row[key]) for key in ("fraction_very_recent_or_extensive", "fraction_recent_or_older", "fraction_no_recent_recombination", "fraction_unresolved")) - 1) < 1e-8 for row in recomb)


def test_recombination_source_concordance_and_classes() -> None:
    root = Path(__file__).resolve().parents[1]
    concordance = read_tsv(root / "results" / "stage2_recombination_source_concordance.tsv")
    assert len(concordance) == 3
    assert all(row["concordance"] == "1" for row in concordance)
    spec = importlib.util.spec_from_file_location("recomb", SCRIPT_DIR / "02l_reconstruct_recombination_track.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    assert module.classify_topology_code(("1", "1", "2")) == "VERY_RECENT_OR_EXTENSIVE"
    assert module.classify_topology_code(("1", "0", "1")) == "RECENT_OR_OLDER"
    assert module.classify_topology_code(("0", "0", "0")) == "NO_RECENT_RECOMBINATION"
    assert module.classify_topology_code(("1", "9", "1")) == "UNRESOLVED"
