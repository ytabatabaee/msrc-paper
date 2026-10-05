#!/usr/bin/env python3
"""Initial descriptive summaries for frozen house-mouse local trees."""

from __future__ import annotations

import argparse
import struct
import sys
import zlib
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from house_mouse_t_complex_utils import (  # noqa: E402
    FIGURES,
    METADATA,
    PROCESSED,
    RESULTS,
    load_mapping,
    parse_newick,
    read_tsv,
    write_tsv,
)


TREE_FILE = PROCESSED / "house_mouse_t_complex_ml_5kb.tre"
TREE_META = PROCESSED / "house_mouse_t_complex_ml_5kb_metadata.tsv"
SUMMARY_TSV = RESULTS / "stage1_local_tree_summary.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    return parser.parse_args()


def split_for_four(tips: list[str], groups: dict[str, str]) -> str:
    target = ["Mus musculus domesticus", "Mus musculus musculus", "Mus musculus castaneus", "Mus spretus"]
    present = {g: [t for t in tips if groups.get(t) == g] for g in target}
    if any(not present[g] for g in target):
        return "NA"
    return "quartet_available_mapping_not_topology_scored"


def write_placeholder_figures(rows: list[dict[str, object]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    x = [float(r["midpoint_bp"]) for r in rows if r.get("midpoint_bp") not in {"NA", None}]
    y = [float(r["fraction_tips_present"]) for r in rows if r.get("midpoint_bp") not in {"NA", None}]
    width, height = 1000, 360
    margin_l, margin_r, margin_t, margin_b = 80, 30, 35, 55
    pixels = bytearray([255] * (width * height * 3))

    def set_pixel(px: int, py: int, color: tuple[int, int, int]) -> None:
        if 0 <= px < width and 0 <= py < height:
            idx = (py * width + px) * 3
            pixels[idx : idx + 3] = bytes(color)

    def line(x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        dx = abs(x1 - x0)
        dy = -abs(y1 - y0)
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        err = dx + dy
        while True:
            for ox in (-1, 0, 1):
                for oy in (-1, 0, 1):
                    set_pixel(x0 + ox, y0 + oy, color)
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                x0 += sx
            if e2 <= dx:
                err += dx
                y0 += sy

    axis = (20, 20, 20)
    line(margin_l, height - margin_b, width - margin_r, height - margin_b, axis)
    line(margin_l, margin_t, margin_l, height - margin_b, axis)
    points: list[tuple[int, int]] = []
    if x:
        xmin, xmax = min(x), max(x)
        span = xmax - xmin if xmax > xmin else 1.0
        for xi, yi in zip(x, y, strict=True):
            px = int(margin_l + (xi - xmin) / span * (width - margin_l - margin_r))
            py = int(height - margin_b - max(0.0, min(1.0, yi)) * (height - margin_t - margin_b))
            points.append((px, py))
        for a, b in zip(points, points[1:]):
            line(a[0], a[1], b[0], b[1], (37, 99, 135))
    raw = b"".join(b"\x00" + pixels[row * width * 3 : (row + 1) * width * 3] for row in range(height))

    def chunk(kind: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)

    png = (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )
    (FIGURES / "house_mouse_t_complex_tree_inventory.png").write_bytes(png)

    pdf_lines = [
        "%PDF-1.4\n",
        "1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj\n",
        "2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj\n",
        "3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 1000 360] /Contents 4 0 R >> endobj\n",
    ]
    commands = ["0.08 0.08 0.08 RG 1 w 80 55 m 970 55 l S 80 55 m 80 325 l S"]
    if points:
        path = [f"{points[0][0]} {height - points[0][1]} m"]
        path.extend(f"{px} {height - py} l" for px, py in points[1:])
        commands.append("0.15 0.39 0.53 RG 1.5 w " + " ".join(path) + " S")
    commands.append("BT /F1 14 Tf 80 335 Td (House mouse t-complex local tree inventory) Tj ET")
    stream = "\n".join(commands).encode("ascii")
    pdf_lines.append(f"4 0 obj << /Length {len(stream)} /Resources << /Font << /F1 5 0 R >> >> >> stream\n")
    head = "".join(pdf_lines).encode("ascii")
    body = stream + b"\nendstream endobj\n5 0 obj << /Type /Font /Subtype /Type1 /BaseFont /Helvetica >> endobj\n"
    offsets = [0]
    parts = [b"%PDF-1.4\n"]
    objects = [
        b"1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj\n",
        b"2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj\n",
        b"3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 1000 360] /Contents 4 0 R /Resources << /Font << /F1 5 0 R >> >> >> endobj\n",
        f"4 0 obj << /Length {len(stream)} >> stream\n".encode("ascii") + stream + b"\nendstream endobj\n",
        b"5 0 obj << /Type /Font /Subtype /Type1 /BaseFont /Helvetica >> endobj\n",
    ]
    parts = [b"%PDF-1.4\n"]
    offsets = []
    for obj in objects:
        offsets.append(sum(len(part) for part in parts))
        parts.append(obj)
    xref_pos = sum(len(part) for part in parts)
    xref = ["xref\n0 6\n0000000000 65535 f \n"]
    xref.extend(f"{offset:010d} 00000 n \n" for offset in offsets)
    trailer = f"trailer << /Size 6 /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n"
    parts.append(("".join(xref) + trailer).encode("ascii"))
    (FIGURES / "house_mouse_t_complex_tree_inventory.pdf").write_bytes(b"".join(parts))


def main() -> int:
    parse_args()
    if not TREE_FILE.exists() or not TREE_META.exists():
        raise SystemExit("Stage 1 gene-tree file and metadata table are required before 01b.")
    lines = [line.strip() for line in TREE_FILE.read_text().splitlines() if line.strip()]
    meta = read_tsv(TREE_META)
    mapping = {row["tree_tip"]: row.get("subspecies", "NA") for row in load_mapping(METADATA / "tip_mapping.tsv")}
    parsed = [parse_newick(line) for line in lines]
    all_tips = sorted({tip for item in parsed for tip in item["tips"]})
    topology_strings = Counter(line for line in lines)
    rows = []
    for row, item, newick in zip(meta, parsed, lines, strict=True):
        tips = set(item["tips"])
        rows.append(
            {
                "locus_id": row["locus_id"],
                "start_bp": row["start_bp"],
                "end_bp": row["end_bp"],
                "midpoint_bp": row["midpoint_bp"],
                "n_tips": len(tips),
                "n_missing_tips": len(set(all_tips) - tips),
                "fraction_tips_present": f"{len(tips) / len(all_tips):.8f}" if all_tips else "NA",
                "topology_exact_newick_count": topology_strings[newick],
                "domesticus_musculus_castaneus_spretus_quartet_status": split_for_four(list(tips), mapping),
                "q1": "NA",
                "q2": "NA",
                "q3": "NA",
            }
        )
    write_tsv(
        SUMMARY_TSV,
        rows,
        [
            "locus_id",
            "start_bp",
            "end_bp",
            "midpoint_bp",
            "n_tips",
            "n_missing_tips",
            "fraction_tips_present",
            "topology_exact_newick_count",
            "domesticus_musculus_castaneus_spretus_quartet_status",
            "q1",
            "q2",
            "q3",
        ],
    )
    write_placeholder_figures(rows)
    print(f"Wrote {SUMMARY_TSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
