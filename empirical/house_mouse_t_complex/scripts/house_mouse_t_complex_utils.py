#!/usr/bin/env python3
"""Shared utilities for the house-mouse t-complex empirical audit."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import subprocess
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
EMPIRICAL_ROOT = REPO_ROOT / "empirical" / "house_mouse_t_complex"
DATA_ROOT = REPO_ROOT / "data" / "house_mouse_t_complex"
RESULTS = EMPIRICAL_ROOT / "results"
FIGURES = EMPIRICAL_ROOT / "figures"
PROCESSED = DATA_ROOT / "processed"
METADATA = DATA_ROOT / "metadata"

PRIMARY_ARCHIVE = (
    "Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/"
    "1-Trees_for_all_5kb_windows/ML_trees.zip"
)
TREE_TOPOLOGY_DIR = (
    "Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/"
    "2-Tree_topologies_for_all_5kb_windows/1-ML_IQtree/"
)
ML_ARCHIVES = {
    "PASS_SNPs": "Data/1-PASS_SNPs/1-Trees_for_all_5kb_windows/ML_trees.zip",
    "Coverage_AlleleRatio_RAW_SNPs": PRIMARY_ARCHIVE,
    "Coverage_RAW_SNPs": "Data/3-Coverage-Filtered_RAW_SNPs/1-Trees_for_all_5kb_windows/ML_trees.zip",
}

TREE_EXTENSIONS = {".treefile", ".contree", ".tre", ".tree", ".nwk", ".newick"}
AUXILIARY_EXTENSIONS = {
    ".bionj",
    ".boottrees",
    ".ckp.gz",
    ".iqtree",
    ".log",
    ".mldist",
    ".model.gz",
    ".splits.nex",
    ".ufboot",
}
TSV_LINETERMINATOR = "\n"


@dataclass(frozen=True)
class TreeRecord:
    window_id: str
    source_filename: str
    newick: str
    sha256: str
    start_bp: int | None
    end_bp: int | None
    midpoint_bp: float | None
    source_start_offset: int | None
    source_end_offset: int | None
    tree_valid: bool
    n_tips: int
    tips: tuple[str, ...]
    duplicate_tips: tuple[str, ...]
    branch_lengths_present: bool
    negative_branch_lengths: bool
    internal_labels_present: bool
    root_representation: str
    parse_error: str


class NewickParseError(ValueError):
    """Raised when a Newick string cannot be parsed conservatively."""


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.i = 0
        self.tips: list[str] = []
        self.internal_labels: list[str] = []
        self.branch_lengths: list[float] = []
        self.branch_lengths_present = False
        self.negative_branch_lengths = False

    def parse(self) -> tuple[list[str], list[str], list[float], int]:
        if not self.text.endswith(";"):
            raise NewickParseError("missing terminal semicolon")
        root_children = self._subtree()
        self._skip_ws()
        if self.i >= len(self.text) or self.text[self.i] != ";":
            raise NewickParseError("expected terminal semicolon")
        self.i += 1
        self._skip_ws()
        if self.i != len(self.text):
            raise NewickParseError("unexpected characters after semicolon")
        return self.tips, self.internal_labels, self.branch_lengths, root_children

    def _subtree(self) -> int:
        self._skip_ws()
        if self.i >= len(self.text):
            raise NewickParseError("unexpected end of string")
        if self.text[self.i] == "(":
            self.i += 1
            n_children = 0
            while True:
                self._subtree()
                n_children += 1
                self._skip_ws()
                if self.i >= len(self.text):
                    raise NewickParseError("unterminated internal node")
                if self.text[self.i] == ",":
                    self.i += 1
                    continue
                if self.text[self.i] == ")":
                    self.i += 1
                    break
                raise NewickParseError(f"expected comma or close parenthesis at byte {self.i}")
            label = self._read_label()
            if label:
                self.internal_labels.append(label)
            self._read_branch_length()
            return n_children
        label = self._read_label()
        if not label:
            raise NewickParseError(f"missing leaf label at byte {self.i}")
        self.tips.append(label)
        self._read_branch_length()
        return 0

    def _skip_ws(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _read_label(self) -> str:
        self._skip_ws()
        if self.i >= len(self.text):
            return ""
        if self.text[self.i] == "'":
            self.i += 1
            chars: list[str] = []
            while self.i < len(self.text):
                ch = self.text[self.i]
                self.i += 1
                if ch == "'":
                    if self.i < len(self.text) and self.text[self.i] == "'":
                        chars.append("'")
                        self.i += 1
                        continue
                    return "".join(chars)
                chars.append(ch)
            raise NewickParseError("unterminated quoted label")
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ":,();":
            self.i += 1
        return self.text[start : self.i].strip()

    def _read_branch_length(self) -> None:
        self._skip_ws()
        if self.i >= len(self.text) or self.text[self.i] != ":":
            return
        self.branch_lengths_present = True
        self.i += 1
        self._skip_ws()
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",();":
            self.i += 1
        token = self.text[start : self.i].strip()
        if not token:
            raise NewickParseError("empty branch length")
        try:
            value = float(token)
        except ValueError as exc:
            raise NewickParseError(f"invalid branch length: {token}") from exc
        self.branch_lengths.append(value)
        if value < 0:
            self.negative_branch_lengths = True


def ensure_dirs() -> None:
    for path in (EMPIRICAL_ROOT, EMPIRICAL_ROOT / "scripts", RESULTS, FIGURES, DATA_ROOT, PROCESSED, METADATA):
        path.mkdir(parents=True, exist_ok=True)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_text_lossy(data: bytes) -> str:
    for encoding in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            pass
    return data.decode("utf-8", errors="replace")


def norm_zip_name(name: str) -> str:
    return name.replace("\\", "/").lstrip("./")


def find_zip_member(zip_file: zipfile.ZipFile, expected: str) -> zipfile.ZipInfo | None:
    expected_norm = norm_zip_name(expected)
    for info in zip_file.infolist():
        name = norm_zip_name(info.filename)
        if name == expected_norm or name.endswith("/" + expected_norm):
            return info
    return None


def nested_zip_bytes(outer: zipfile.ZipFile, expected: str) -> tuple[bytes | None, zipfile.ZipInfo | None]:
    info = find_zip_member(outer, expected)
    if info is None:
        return None, None
    return outer.read(info), info


def extension_for_zip_name(name: str) -> str:
    lower = norm_zip_name(name).lower()
    for suffix in sorted(AUXILIARY_EXTENSIONS | TREE_EXTENSIONS, key=len, reverse=True):
        if lower.endswith(suffix):
            return suffix
    return Path(lower).suffix


def is_tree_candidate(name: str) -> bool:
    return extension_for_zip_name(name) in TREE_EXTENSIONS


def final_tree_files(names: Iterable[str]) -> list[str]:
    files = sorted(n for n in names if not n.endswith("/") and is_tree_candidate(n))
    treefiles = [n for n in files if extension_for_zip_name(n) == ".treefile"]
    if treefiles:
        return treefiles
    return files


def parse_source_offsets_from_name(name: str) -> tuple[int | None, int | None, str]:
    base = Path(norm_zip_name(name)).name
    dash_match = re.search(r"(?<!\d)(\d+)-(\d+)(?=\D)", base)
    if dash_match:
        start, end = int(dash_match.group(1)), int(dash_match.group(2))
        if end < start:
            return None, None, "unusable"
        return start, end, "explicit_half_open_start_end"
    integers = [int(x) for x in re.findall(r"(?<![A-Za-z])(\d{4,})(?![A-Za-z])", base)]
    if len(integers) >= 2:
        start, end = integers[0], integers[1]
        if end < start:
            return None, None, "unusable"
        return start, end, "explicit_start_end_unverified_interval_type"
    if len(integers) == 1 and re.search(r"(?:start|pos|bp|window|chr17|5kb)", base, re.IGNORECASE):
        start = integers[0]
        return start, start + 5000, "single_start_plus_5000_half_open"
    return None, None, "unresolved"


def parse_coordinates_from_name(name: str, chr17_offset: int = 5_000_000) -> tuple[int | None, int | None, str]:
    start_offset, end_offset, status = parse_source_offsets_from_name(name)
    if start_offset is None or end_offset is None:
        return None, None, status
    if end_offset - start_offset == 5000:
        return chr17_offset + start_offset, chr17_offset + end_offset - 1, status + "_chr17_plus_5000000"
    return chr17_offset + start_offset, chr17_offset + end_offset, status + "_chr17_plus_5000000"


def canonical_tip_label(raw_tip: str, source_start_offset: int | None, source_end_offset: int | None) -> str:
    if source_start_offset is None or source_end_offset is None:
        return raw_tip
    suffix = f"_{source_start_offset}-{source_end_offset}"
    if raw_tip.endswith(suffix):
        return raw_tip[: -len(suffix)]
    return raw_tip


def canonicalize_newick_labels(text: str, raw_tips: Iterable[str], source_start_offset: int | None, source_end_offset: int | None) -> str:
    canonical = text
    replacements = {
        raw: canonical_tip_label(raw, source_start_offset, source_end_offset)
        for raw in raw_tips
        if raw != canonical_tip_label(raw, source_start_offset, source_end_offset)
    }
    for raw in sorted(replacements, key=len, reverse=True):
        canonical = canonical.replace(raw, replacements[raw])
    return canonical


def parse_newick(text: str) -> dict[str, object]:
    parser = NewickParser(text)
    tips, internal, branches, root_children = parser.parse()
    duplicates = tuple(sorted(label for label, n in Counter(tips).items() if n > 1))
    if root_children == 2:
        root_representation = "rooted_or_binary_root"
    elif root_children > 2:
        root_representation = "unrooted_or_polytomy_root"
    else:
        root_representation = "single_tip_or_unknown"
    return {
        "tips": tuple(tips),
        "duplicate_tips": duplicates,
        "branch_lengths_present": parser.branch_lengths_present,
        "negative_branch_lengths": parser.negative_branch_lengths,
        "internal_labels_present": bool(internal),
        "internal_labels": tuple(internal),
        "root_representation": root_representation,
        "n_branch_lengths": len(branches),
    }


def make_tree_records(nested_data: bytes) -> tuple[list[TreeRecord], dict[str, object]]:
    records: list[TreeRecord] = []
    with zipfile.ZipFile(io.BytesIO(nested_data)) as nested:
        names = sorted(info.filename for info in nested.infolist() if not info.is_dir())
        finals = final_tree_files(names)
        for index, name in enumerate(finals, start=1):
            raw = nested.read(name)
            text = read_text_lossy(raw).strip()
            source_start, source_end, _source_coord_status = parse_source_offsets_from_name(name)
            start, end, coord_status = parse_coordinates_from_name(name)
            midpoint = (start + end) / 2 if start is not None and end is not None else None
            valid = True
            parse_error = ""
            parsed: dict[str, object] = {}
            try:
                parsed = parse_newick(text)
                raw_tips = tuple(parsed["tips"])
                canonical_tips = tuple(canonical_tip_label(t, source_start, source_end) for t in raw_tips)
                duplicate_canonical = tuple(sorted(label for label, n in Counter(canonical_tips).items() if n > 1))
                parsed["raw_tips"] = raw_tips
                parsed["tips"] = canonical_tips
                parsed["duplicate_tips"] = duplicate_canonical
                text = canonicalize_newick_labels(text, raw_tips, source_start, source_end)
            except Exception as exc:  # noqa: BLE001 - audit records exact failure, then continues.
                valid = False
                parse_error = str(exc)
                parsed = {
                    "tips": tuple(),
                    "duplicate_tips": tuple(),
                    "branch_lengths_present": False,
                    "negative_branch_lengths": False,
                    "internal_labels_present": False,
                    "root_representation": "malformed",
                }
            records.append(
                TreeRecord(
                    window_id=f"w{index:06d}",
                    source_filename=name,
                    newick=text,
                    sha256=sha256_bytes(raw),
                    start_bp=start,
                    end_bp=end,
                    midpoint_bp=midpoint,
                    source_start_offset=source_start,
                    source_end_offset=source_end,
                    tree_valid=valid,
                    n_tips=len(parsed["tips"]),
                    tips=tuple(parsed["tips"]),
                    duplicate_tips=tuple(parsed["duplicate_tips"]),
                    branch_lengths_present=bool(parsed["branch_lengths_present"]),
                    negative_branch_lengths=bool(parsed["negative_branch_lengths"]),
                    internal_labels_present=bool(parsed["internal_labels_present"]),
                    root_representation=str(parsed["root_representation"]),
                    parse_error=parse_error,
                )
            )
    summary = coordinate_summary(records)
    summary["coordinate_recovery_status"] = coord_status if records else "no_tree_files"
    return records, summary


def coordinate_summary(records: list[TreeRecord]) -> dict[str, object]:
    valid_coords = [r for r in records if r.start_bp is not None and r.end_bp is not None]
    starts = [r.start_bp for r in valid_coords if r.start_bp is not None]
    coord_recoverable = len(valid_coords) == len(records) and bool(records)
    duplicates = sorted(k for k, v in Counter((r.start_bp, r.end_bp) for r in valid_coords).items() if v > 1)
    missing: list[str] = []
    non_overlapping = None
    observed_step = None
    if coord_recoverable:
        ordered = sorted(valid_coords, key=lambda r: (r.start_bp or -1, r.end_bp or -1, r.source_filename))
        spans = [(r.start_bp, r.end_bp) for r in ordered]
        steps = [(ordered[i + 1].start_bp or 0) - (ordered[i].start_bp or 0) for i in range(len(ordered) - 1)]
        widths = [(r.end_bp or 0) - (r.start_bp or 0) + 1 for r in ordered]
        observed_step = sorted(Counter(steps).items(), key=lambda x: (-x[1], x[0]))[0][0] if steps else None
        non_overlapping = all(spans[i][1] < spans[i + 1][0] for i in range(len(spans) - 1))
        unique_starts = sorted({r.start_bp for r in ordered if r.start_bp is not None})
        if widths and all(width == 5000 for width in widths) and unique_starts:
            expected = set(range(unique_starts[0], unique_starts[-1] + 1, 5000))
            observed = set(unique_starts)
            missing = [str(x) for x in sorted(expected - observed)]
    return {
        "coordinates_recoverable_for_all_final_trees": coord_recoverable,
        "n_coordinate_resolved": len(valid_coords),
        "observed_start_bp_min": min(starts) if starts else None,
        "observed_start_bp_max": max(starts) if starts else None,
        "observed_window_count": len(records),
        "observed_modal_step_bp": observed_step,
        "non_overlapping_by_coordinates": non_overlapping,
        "duplicate_windows": ["-".join(map(str, pair)) for pair in duplicates],
        "missing_window_starts": missing,
        "coordinate_basis": "chr17 coordinate = archive window coordinate + 5,000,000, from Data/Readme.txt",
    }


def tree_inventory(records: list[TreeRecord]) -> dict[str, object]:
    valid = [r for r in records if r.tree_valid]
    malformed = [r for r in records if not r.tree_valid]
    tip_counter: Counter[str] = Counter()
    missing_counter: Counter[str] = Counter()
    duplicate_files = [r.source_filename for r in records if r.duplicate_tips]
    for r in valid:
        tip_counter.update(set(r.tips))
    union = sorted(tip_counter)
    intersection = sorted(t for t in union if tip_counter[t] == len(valid))
    for tip in union:
        missing_counter[tip] = len(valid) - tip_counter[tip]
    return {
        "n_final_tree_files": len(records),
        "n_valid_newicks": len(valid),
        "n_malformed_newicks": len(malformed),
        "tip_count_distribution": dict(sorted(Counter(r.n_tips for r in valid).items())),
        "root_representation_distribution": dict(sorted(Counter(r.root_representation for r in valid).items())),
        "union_tip_labels": union,
        "intersection_tip_labels": intersection,
        "missing_tip_frequency": dict(sorted(missing_counter.items())),
        "files_with_duplicate_tip_labels": duplicate_files,
        "branch_lengths_present_all_valid": all(r.branch_lengths_present for r in valid) if valid else False,
        "branch_lengths_present_any_valid": any(r.branch_lengths_present for r in valid),
        "support_or_internal_labels_present_any_valid": any(r.internal_labels_present for r in valid),
        "negative_branch_lengths_any_valid": any(r.negative_branch_lengths for r in valid),
        "malformed_files": [{"source_filename": r.source_filename, "parse_error": r.parse_error} for r in malformed],
    }


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator=TSV_LINETERMINATOR)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_tsv_value(row.get(field)) for field in fields})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def format_tsv_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, (list, tuple, set)):
        return ",".join(str(x) for x in value) if value else ""
    return str(value)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:  # noqa: BLE001
        return "unavailable"


def tip_inventory_rows(records: list[TreeRecord]) -> list[dict[str, object]]:
    valid = [r for r in records if r.tree_valid]
    n_valid = len(valid)
    counts: Counter[str] = Counter()
    for r in valid:
        counts.update(set(r.tips))
    rows = []
    for tip in sorted(counts):
        rows.append(
            {
                "tip_label": tip,
                "n_trees_present": counts[tip],
                "fraction_trees_present": f"{counts[tip] / n_valid:.8f}" if n_valid else "NA",
                "inferred_category": "unknown",
                "evidence": "No source-supported tip-to-subspecies/t-status mapping applied by Stage 0.",
                "mapping_status": "unresolved",
            }
        )
    return rows


def window_inventory_rows(records: list[TreeRecord]) -> list[dict[str, object]]:
    ordered = sorted(records, key=lambda r: ((r.start_bp is None, r.start_bp or 10**30), r.source_filename))
    return [
        {
            "window_id": r.window_id,
            "source_filename": r.source_filename,
            "start_bp": r.start_bp,
            "end_bp": r.end_bp,
            "midpoint_bp": r.midpoint_bp,
            "source_start_offset": r.source_start_offset,
            "source_end_offset": r.source_end_offset,
            "tree_valid": r.tree_valid,
            "n_tips": r.n_tips,
            "sha256": r.sha256,
        }
        for r in ordered
    ]


def default_tip_mapping_rows(tips: Iterable[str]) -> list[dict[str, object]]:
    return [infer_tip_mapping(tip) for tip in sorted(tips)]


def infer_tip_mapping(tip: str) -> dict[str, object]:
    evidence = (
        "Kelemen & Vicoso 2018 Genetics article states AFG/CZE/KAZ are M. m. musculus, "
        "GER/FRA are M. m. domesticus, CAST is M. m. castaneus, and SPRE is M. spretus; "
        "Data/Readme.txt and tree labels identify t-haplotype subsets."
    )
    row = {
        "tree_tip": tip,
        "individual": tip,
        "subspecies": "NA",
        "population": "NA",
        "t_status": "NA",
        "mapping_source": evidence,
        "mapping_confidence": "unresolved",
        "notes": "",
    }
    rules = [
        ("Mmc_CAST", "Mus musculus castaneus", "CAST"),
        ("Mmd_FRA", "Mus musculus domesticus", "France"),
        ("Mmd_GER", "Mus musculus domesticus", "Germany"),
        ("Mmm_AFG", "Mus musculus musculus", "Afghanistan"),
        ("Mmm_CZE", "Mus musculus musculus", "Czech Republic"),
        ("Mmm_KAZ", "Mus musculus musculus", "Kazakhstan"),
        ("Ms_SPRE", "Mus spretus", "SPRE"),
    ]
    for prefix, subspecies, population in rules:
        if tip.startswith(prefix):
            row["subspecies"] = subspecies
            row["population"] = population
            row["mapping_confidence"] = "strong"
            break
    if row["mapping_confidence"] == "unresolved":
        row["notes"] = "No source-backed prefix rule matched this label."
        return row
    if "_tHaplSubset.fa" in tip:
        row["t_status"] = "pseudo-t_haplotype"
        row["notes"] = "Tree tip is the pseudo-t/t-haplotype subset sequence for a heterozygous t-carrier."
    elif str(row["subspecies"]) == "Mus spretus":
        row["t_status"] = "outgroup_not_t_haplotype"
        row["notes"] = "M. spretus is the outgroup; the source paper states t-haplotypes are in the M. musculus species complex, not M. spretus."
    else:
        row["t_status"] = "standard_noncarrier"
        row["notes"] = "Tree tip is a noncarrier/standard sequence under the source paper's 15 t-haplotypes plus 40 noncarrier-tree design."
    return row


def load_mapping(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    return read_tsv(path)


def mapping_gate(mapping_rows: list[dict[str, str]], required_tips: Iterable[str]) -> dict[str, object]:
    by_tip = {row.get("tree_tip", ""): row for row in mapping_rows}
    missing_rows = sorted(set(required_tips) - set(by_tip))
    usable = []
    unresolved = []
    tentative = []
    for tip in sorted(required_tips):
        row = by_tip.get(tip)
        if not row:
            continue
        confidence = row.get("mapping_confidence", "unresolved")
        subspecies = row.get("subspecies", "NA")
        if confidence == "tentative":
            tentative.append(tip)
        if confidence in {"exact", "strong"} and subspecies and subspecies != "NA":
            usable.append(tip)
        else:
            unresolved.append(tip)
    groups = sorted({by_tip[t].get("subspecies", "") for t in usable if by_tip[t].get("subspecies", "") not in {"", "NA"}})
    t_status_values = sorted(
        {
            by_tip[t].get("t_status", "")
            for t in usable
            if by_tip[t].get("t_status", "") not in {"", "NA", "unresolved"}
        }
    )
    enough_groups = len(groups) >= 4
    mostly_resolved = len(usable) > 0 and len(usable) / max(1, len(set(required_tips))) >= 0.8 and enough_groups
    t_status_resolved = len(t_status_values) >= 2 and len(usable) == len(set(required_tips))
    return {
        "missing_mapping_rows": missing_rows,
        "usable_exact_or_strong_tips": sorted(usable),
        "unresolved_or_tentative_tips": sorted(set(unresolved + tentative)),
        "usable_subspecies_groups": groups,
        "t_status_values": t_status_values,
        "t_status_resolved": t_status_resolved,
        "n_usable_exact_or_strong_tips": len(usable),
        "n_required_tips": len(set(required_tips)),
        "mostly_resolved": mostly_resolved,
        "astral_ready": mostly_resolved,
    }


def metadata_candidate_names(names: Iterable[str]) -> list[str]:
    patterns = re.compile(
        r"(readme|metadata|sample|legend|supplement|supplementary|table|code|script|mapping|individual|population)",
        re.IGNORECASE,
    )
    candidates = []
    for name in sorted(names):
        lower = name.lower()
        if name.endswith("/"):
            continue
        if patterns.search(name) or lower.endswith((".txt", ".tsv", ".csv", ".xls", ".xlsx", ".pdf", ".r", ".sh", ".py")):
            candidates.append(name)
    return candidates


def topology_result_summary(outer: zipfile.ZipFile) -> dict[str, object]:
    rows = []
    prefix = TREE_TOPOLOGY_DIR
    for info in sorted(outer.infolist(), key=lambda x: x.filename):
        name = norm_zip_name(info.filename)
        if prefix not in name or info.is_dir():
            continue
        if "Tree_results_" not in Path(name).name:
            continue
        sample = read_text_lossy(outer.read(info))[:2000]
        header = sample.splitlines()[0] if sample.splitlines() else ""
        rows.append(
            {
                "filename": name,
                "size_bytes": info.file_size,
                "first_line": header,
                "contains_tip_like_labels": bool(re.search(r"[A-Za-z]+[_-][A-Za-z0-9]+", sample)),
                "contains_coordinates": bool(re.search(r"\d{5,}", sample)),
            }
        )
    return {"n_tree_result_files": len(rows), "files": rows[:50]}


def archive_audit(archive: Path, *, write_outputs: bool = True) -> dict[str, object]:
    if write_outputs:
        ensure_dirs()
    archive = archive.resolve()
    source_size = archive.stat().st_size
    source_sha = sha256_path(archive)
    with zipfile.ZipFile(archive) as outer:
        outer_names = sorted(info.filename for info in outer.infolist())
        nested_status = {}
        for label, expected in ML_ARCHIVES.items():
            info = find_zip_member(outer, expected)
            nested_status[label] = {
                "expected_path": expected,
                "present": info is not None,
                "size_bytes": info.file_size if info else None,
                "archive_member": norm_zip_name(info.filename) if info else None,
            }
        primary_data, primary_info = nested_zip_bytes(outer, PRIMARY_ARCHIVE)
        metadata_candidates = metadata_candidate_names(outer_names)
        topology_summary = topology_result_summary(outer)
    records: list[TreeRecord] = []
    nested_summary: dict[str, object] = {"available": False}
    coord_summary: dict[str, object] = {}
    inv: dict[str, object] = {}
    nested_zip_sha = None
    if primary_data is not None:
        nested_zip_sha = sha256_bytes(primary_data)
        with zipfile.ZipFile(io.BytesIO(primary_data)) as nested:
            nested_names = sorted(info.filename for info in nested.infolist() if not info.is_dir())
            final_names = final_tree_files(nested_names)
            ext_counts = dict(sorted(Counter(extension_for_zip_name(n) for n in nested_names).items()))
            pattern_counts = {
                "final_tree_candidates": len(final_names),
                "treefile": sum(extension_for_zip_name(n) == ".treefile" for n in nested_names),
                "auxiliary_or_other": len(nested_names) - len(final_names),
            }
            nested_summary = {
                "available": True,
                "primary_member": norm_zip_name(primary_info.filename) if primary_info else PRIMARY_ARCHIVE,
                "primary_nested_zip_sha256": nested_zip_sha,
                "number_of_files": len(nested_names),
                "extensions": ext_counts,
                "filename_patterns": pattern_counts,
                "has_multiple_iqtree_auxiliary_files_per_window": bool(
                    pattern_counts["auxiliary_or_other"] and pattern_counts["treefile"]
                ),
                "first_20_representative_filenames": nested_names[:20],
                "last_20_representative_filenames": nested_names[-20:],
            }
        records, coord_summary = make_tree_records(primary_data)
        inv = tree_inventory(records)
    tip_rows = tip_inventory_rows(records)
    window_rows = window_inventory_rows(records)
    mapping_path = METADATA / "tip_mapping.tsv"
    if write_outputs:
        write_tsv(
            PROCESSED / "stage0_tree_tip_inventory.tsv",
            tip_rows,
        ["tip_label", "n_trees_present", "fraction_trees_present", "inferred_category", "evidence", "mapping_status"],
        )
        write_tsv(
            PROCESSED / "stage0_window_inventory.tsv",
            window_rows,
        ["window_id", "source_filename", "start_bp", "end_bp", "midpoint_bp", "tree_valid", "n_tips", "sha256"],
        )
    existing_mapping_rows = load_mapping(mapping_path) if mapping_path.exists() else []
    generated_unresolved = bool(existing_mapping_rows) and all(
        row.get("mapping_confidence") == "unresolved"
        and row.get("mapping_source") in {"unresolved_from_archive_audit", ""}
        for row in existing_mapping_rows
    )
    if write_outputs and (not mapping_path.exists() or not existing_mapping_rows or generated_unresolved):
        write_tsv(
            mapping_path,
            default_tip_mapping_rows(inv.get("union_tip_labels", [])),
            [
                "tree_tip",
                "individual",
                "subspecies",
                "population",
                "t_status",
                "mapping_source",
                "mapping_confidence",
                "notes",
            ],
        )
    mapping_rows = load_mapping(mapping_path) if write_outputs else default_tip_mapping_rows(inv.get("union_tip_labels", []))
    gate = mapping_gate(mapping_rows, inv.get("union_tip_labels", []))
    stage0_status = stage0_gate_status(primary_data is not None, coord_summary, gate)
    manifest = {
        "analysis": "house_mouse_t_complex_stage0_source_audit",
        "command": " ".join(["00_audit_house_mouse_data.py", "--archive", str(archive)]),
        "git_commit": git_commit(),
        "source_archive": str(archive),
        "source_archive_size_bytes": source_size,
        "source_archive_sha256": source_sha,
        "nested_ml_archives": nested_status,
        "primary_nested_zip_sha256": nested_zip_sha,
        "primary_archive_summary": nested_summary,
        "coordinate_summary": coord_summary,
        "tree_inventory": inv,
        "topology_result_summary": topology_summary,
        "metadata_candidates": metadata_candidates[:200],
        "mapping_gate": gate,
        "stage0_gate": stage0_status,
        "outputs": {
            "tip_inventory": str((PROCESSED / "stage0_tree_tip_inventory.tsv").relative_to(REPO_ROOT)),
            "window_inventory": str((PROCESSED / "stage0_window_inventory.tsv").relative_to(REPO_ROOT)),
            "tip_mapping": str(mapping_path.relative_to(REPO_ROOT)),
        },
    }
    if write_outputs:
        write_json(RESULTS / "stage0_manifest.json", manifest)
        (RESULTS / "stage0_source_audit.md").write_text(stage0_markdown(manifest) + "\n")
    return manifest


def stage0_gate_status(primary_available: bool, coord_summary: dict[str, object], gate: dict[str, object]) -> dict[str, object]:
    a = primary_available
    b = bool(coord_summary.get("coordinates_recoverable_for_all_final_trees"))
    c = bool(gate.get("mostly_resolved"))
    d = bool(gate.get("t_status_resolved"))
    e = bool(gate.get("astral_ready"))
    f = not a
    if not a or not c:
        overall = "FAIL"
    elif not b or not d:
        overall = "PASS_WITH_LIMITATIONS"
    else:
        overall = "PASS"
    return {
        "A_per_window_ML_Newick_trees_available": a,
        "B_windows_mapped_to_coordinates": b,
        "C_tips_mapped_unambiguously_to_subspecies": c,
        "D_t_status_distinguishable": d,
        "E_enough_for_ASTRAL_mapping": e,
        "F_obvious_blocker": f,
        "overall_status": overall,
    }


def stage0_markdown(manifest: dict[str, object]) -> str:
    gate = manifest["stage0_gate"]
    inv = manifest.get("tree_inventory", {})
    coord = manifest.get("coordinate_summary", {})
    nested = manifest.get("nested_ml_archives", {})
    primary = manifest.get("primary_archive_summary", {})
    lines = [
        "# House mouse t-complex Stage 0 source audit",
        "",
        "This audit inspects the Kelemen & Vicoso archive without requiring full manual extraction.",
        "",
        "## Source archive",
        "",
        f"- Archive: `{manifest['source_archive']}`",
        f"- Size bytes: {manifest['source_archive_size_bytes']}",
        f"- SHA256: `{manifest['source_archive_sha256']}`",
        f"- Git commit: `{manifest['git_commit']}`",
        "",
        "## Nested ML tree archives",
        "",
    ]
    for label, row in nested.items():
        lines.append(
            f"- {label}: present={row['present']}; size_bytes={row['size_bytes']}; member=`{row['archive_member']}`"
        )
    lines.extend(
        [
            "",
            "## Primary nested archive",
            "",
            f"- Available: {primary.get('available')}",
            f"- Number of files: {primary.get('number_of_files', 'NA')}",
            f"- Extensions: `{json.dumps(primary.get('extensions', {}), sort_keys=True)}`",
            f"- Filename patterns: `{json.dumps(primary.get('filename_patterns', {}), sort_keys=True)}`",
            f"- Multiple IQ-TREE auxiliary files per window: {primary.get('has_multiple_iqtree_auxiliary_files_per_window', 'NA')}",
            f"- First 20 filenames: `{json.dumps(primary.get('first_20_representative_filenames', []))}`",
            f"- Last 20 filenames: `{json.dumps(primary.get('last_20_representative_filenames', []))}`",
            "",
            "## Coordinate audit",
            "",
            f"- Coordinates recoverable for all final trees: {coord.get('coordinates_recoverable_for_all_final_trees', False)}",
            f"- Resolved windows: {coord.get('n_coordinate_resolved', 0)} / {coord.get('observed_window_count', 0)}",
            f"- Observed start range: {coord.get('observed_start_bp_min', 'NA')} - {coord.get('observed_start_bp_max', 'NA')}",
            f"- Modal start step: {coord.get('observed_modal_step_bp', 'NA')}",
            f"- Non-overlapping by coordinates: {coord.get('non_overlapping_by_coordinates', 'NA')}",
            f"- Duplicate windows: `{json.dumps(coord.get('duplicate_windows', []))}`",
            f"- Missing window starts: `{json.dumps(coord.get('missing_window_starts', []))}`",
            "",
            "## Newick/tree audit",
            "",
            f"- Final tree files: {inv.get('n_final_tree_files', 0)}",
            f"- Valid Newicks: {inv.get('n_valid_newicks', 0)}",
            f"- Malformed Newicks: {inv.get('n_malformed_newicks', 0)}",
            f"- Tip-count distribution: `{json.dumps(inv.get('tip_count_distribution', {}), sort_keys=True)}`",
            f"- Root representation distribution: `{json.dumps(inv.get('root_representation_distribution', {}), sort_keys=True)}`",
            f"- Union tip labels: `{json.dumps(inv.get('union_tip_labels', []))}`",
            f"- Intersection tip labels: `{json.dumps(inv.get('intersection_tip_labels', []))}`",
            f"- Duplicate tip-label files: `{json.dumps(inv.get('files_with_duplicate_tip_labels', []))}`",
            f"- Branch lengths present in all valid trees: {inv.get('branch_lengths_present_all_valid', False)}",
            f"- Internal labels/support values present in any valid tree: {inv.get('support_or_internal_labels_present_any_valid', False)}",
            f"- Negative branch lengths present in any valid tree: {inv.get('negative_branch_lengths_any_valid', False)}",
            "",
            "## Topology result and metadata audit",
            "",
            "The machine-readable manifest records candidate README, metadata, sample, code, legend, and supplementary files found in the archive. Tree_results files are summarized without using them to infer identities unless the source text directly supports a mapping.",
            "",
            f"- Tree_results files found: {manifest.get('topology_result_summary', {}).get('n_tree_result_files', 0)}",
            f"- Metadata candidates recorded: {len(manifest.get('metadata_candidates', []))}",
            "",
            "## Stage 0 scientific gate",
            "",
            f"A. Are the per-window ML Newick trees actually available? **{gate['A_per_window_ML_Newick_trees_available']}**",
            f"B. Can tree windows be mapped to genomic coordinates? **{gate['B_windows_mapped_to_coordinates']}**",
            f"C. Can tips be mapped unambiguously to subspecies? **{gate['C_tips_mapped_unambiguously_to_subspecies']}**",
            f"D. Can t-haplotype/pseudo-t tips be distinguished from standard/noncarrier tips? **{gate['D_t_status_distinguishable']}**",
            f"E. Is there enough information to construct ASTRAL multi-individual/species mappings? **{gate['E_enough_for_ASTRAL_mapping']}**",
            f"F. Is there any obvious reason the dataset cannot support the proposed MSRC analysis? **{gate['F_obvious_blocker']}**",
            "",
            f"Overall status: **{gate['overall_status']}**",
            "",
            "Stage 1 must not be interpreted biologically unless the tip-mapping gate is resolved from source-backed metadata. Tentative mappings are not used automatically for ASTRAL.",
        ]
    )
    return "\n".join(lines)


def stage1_can_run(manifest_path: Path = RESULTS / "stage0_manifest.json") -> tuple[bool, str, dict[str, object] | None]:
    if not manifest_path.exists():
        return False, "Stage 0 manifest is missing.", None
    manifest = json.loads(manifest_path.read_text())
    gate = manifest.get("stage0_gate", {})
    if gate.get("overall_status") not in {"PASS", "PASS_WITH_LIMITATIONS"}:
        return False, f"Stage 0 gate is {gate.get('overall_status')}; Stage 1 is blocked.", manifest
    if not gate.get("C_tips_mapped_unambiguously_to_subspecies"):
        return False, "Stage 0 tip-to-subspecies mapping gate did not pass.", manifest
    return True, "Stage 1 gate passed.", manifest
