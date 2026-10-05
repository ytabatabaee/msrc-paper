"""Lightweight Stage-2 tests that do not require the upstream archive."""

from __future__ import annotations

import sys
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
