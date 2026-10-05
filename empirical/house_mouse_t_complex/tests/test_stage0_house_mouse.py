from __future__ import annotations

import importlib.util
import sys
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
UTILS = ROOT / "empirical" / "house_mouse_t_complex" / "scripts" / "house_mouse_t_complex_utils.py"


spec = importlib.util.spec_from_file_location("house_mouse_t_complex_utils", UTILS)
utils = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules["house_mouse_t_complex_utils"] = utils
spec.loader.exec_module(utils)


def test_nested_zip_extraction_and_tip_inventory(tmp_path: Path) -> None:
    archive = tmp_path / "IST-2017-78-v1+1_Data.zip"
    nested = tmp_path / "ML_trees.zip"
    with zipfile.ZipFile(nested, "w") as z:
        z.writestr("0-5000.fa.contree", "(a_0-5000:1,b_0-5000:1,c_0-5000:1,d_0-5000:1);")
        z.writestr("5000-10000.fa.contree", "(a_5000-10000:1,b_5000-10000:1,c_5000-10000:1);")
    with zipfile.ZipFile(archive, "w") as z:
        z.write(nested, utils.PRIMARY_ARCHIVE)
    manifest = utils.archive_audit(archive, write_outputs=False)
    assert manifest["tree_inventory"]["n_final_tree_files"] == 2
    assert manifest["tree_inventory"]["missing_tip_frequency"]["d"] == 1
    assert manifest["coordinate_summary"]["coordinates_recoverable_for_all_final_trees"] is True


def test_duplicate_missing_window_detection(tmp_path: Path) -> None:
    nested = tmp_path / "ML_trees.zip"
    with zipfile.ZipFile(nested, "w") as z:
        z.writestr("0-5000.fa.contree", "(a_0-5000:1,b_0-5000:1,c_0-5000:1,d_0-5000:1);")
        z.writestr("copy_0-5000.fa.contree", "(a_0-5000:1,b_0-5000:1,c_0-5000:1,d_0-5000:1);")
        z.writestr("10000-15000.fa.contree", "(a_10000-15000:1,b_10000-15000:1,c_10000-15000:1,d_10000-15000:1);")
    records, coord = utils.make_tree_records(nested.read_bytes())
    assert coord["duplicate_windows"] == ["5000000-5004999"]
    assert "5005000" in coord["missing_window_starts"]
    assert len(records) == 3


def test_malformed_newick_and_duplicate_tips(tmp_path: Path) -> None:
    fixture_path = tmp_path / "house_mouse_malformed_fixture.zip"
    with zipfile.ZipFile(fixture_path, "w") as z:
        z.writestr("chr17_5000000_5004999.treefile", "(a:1,a:1,c:1,d:1);")
        z.writestr("chr17_5005000_5009999.treefile", "(a:1,b:1,c:1,d:1)")
    fixture = fixture_path.read_bytes()
    records, _coord = utils.make_tree_records(fixture)
    inv = utils.tree_inventory(records)
    assert records[0].duplicate_tips == ("a",)
    assert inv["n_valid_newicks"] == 1
    assert inv["n_malformed_newicks"] == 1


def test_mapping_confidence_gate() -> None:
    rows = [
        {"tree_tip": "a", "subspecies": "Mus musculus domesticus", "mapping_confidence": "exact"},
        {"tree_tip": "b", "subspecies": "Mus musculus musculus", "mapping_confidence": "strong"},
        {"tree_tip": "c", "subspecies": "Mus musculus castaneus", "mapping_confidence": "tentative"},
        {"tree_tip": "d", "subspecies": "Mus spretus", "mapping_confidence": "unresolved"},
    ]
    gate = utils.mapping_gate(rows, {"a", "b", "c", "d"})
    assert gate["mostly_resolved"] is False
    assert gate["usable_subspecies_groups"] == ["Mus musculus domesticus", "Mus musculus musculus"]
