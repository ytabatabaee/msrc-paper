#!/usr/bin/env python3
"""Build the frozen cross-dataset quartet-mixture/CU distortion result."""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "empirical/cross_dataset_cu"
RESULTS = OUT / "results"
FIGURES = OUT / "figures"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def canonical_split(left: set[str], all_taxa: set[str]) -> str:
    right = all_taxa - left
    a = tuple(sorted(left))
    b = tuple(sorted(right))
    chosen = a if (len(a), a) <= (len(b), b) else b
    other = b if chosen == a else a
    return ",".join(chosen) + "|" + ",".join(other)


def parse_split(text: str) -> str:
    left, right = text.split("|", 1)
    return canonical_split(set(left.split(",")), set(left.split(",")) | set(right.split(",")))


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip().rstrip(";")
        self.i = 0
        self.records: list[tuple[set[str], dict[str, float]]] = []

    def parse(self) -> set[str]:
        leaves = self._subtree()
        self._skip_space()
        if self.i != len(self.text):
            raise ValueError(f"unparsed Newick suffix at {self.i}: {self.text[self.i:self.i+40]}")
        return leaves

    def _skip_space(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _subtree(self) -> set[str]:
        self._skip_space()
        if self.text[self.i] == "(":
            self.i += 1
            descendants: set[str] = set()
            while True:
                descendants |= self._subtree()
                self._skip_space()
                if self.text[self.i] == ",":
                    self.i += 1
                    continue
                if self.text[self.i] == ")":
                    self.i += 1
                    break
                raise ValueError(f"expected comma or close at {self.i}")
            self._skip_space()
            annotation: dict[str, float] = {}
            if self.i < len(self.text) and self.text[self.i] == "'":
                self.i += 1
                start = self.i
                end = self.text.index("'", start)
                annotation = {
                    key: float(value.strip("[]"))
                    for key, value in re.findall(r"([A-Za-z][A-Za-z0-9_]*)=([^;]+)", self.text[start:end])
                }
                self.i = end + 1
            self._skip_length()
            if annotation:
                self.records.append((descendants, annotation))
            return descendants
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",():;'":
            self.i += 1
        label = self.text[start:self.i].strip()
        if not label:
            raise ValueError(f"empty leaf at {self.i}")
        self._skip_length()
        return {label}

    def _skip_length(self) -> None:
        self._skip_space()
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            while self.i < len(self.text) and self.text[self.i] not in ",()":
                self.i += 1


def q1_by_split(path: Path) -> dict[str, float]:
    parser = NewickParser(path.read_text())
    all_taxa = parser.parse()
    out: dict[str, float] = {}
    for leaves, annotation in parser.records:
        if "q1" not in annotation:
            continue
        key = canonical_split(leaves, all_taxa)
        if key in out:
            raise ValueError(f"duplicate annotated split in {path}: {key}")
        out[key] = annotation["q1"]
    return out


def msc_length(q: float) -> float:
    if q <= 1 / 3:
        raise ValueError(f"MSC transform undefined for q <= 1/3: {q}")
    return -math.log(1.5 * (1 - q))


def row(dataset: str, branch_id: str, display: str, n_bg: int, n_aff: int, q_bg: float,
        q_aff: float, q_obs: float, cu_bg: float, cu_all: float) -> dict[str, object]:
    eps = n_aff / (n_bg + n_aff)
    q_pred = (1 - eps) * q_bg + eps * q_aff
    cu_bg_pred = msc_length(q_bg)
    cu_all_pred = msc_length(q_pred)
    return {
        "dataset": dataset,
        "branch_id": branch_id,
        "display_label": display,
        "n_background": n_bg,
        "n_affected": n_aff,
        "epsilon": eps,
        "q_background": q_bg,
        "q_affected": q_aff,
        "q_all_predicted": q_pred,
        "q_all_observed": q_obs,
        "q_prediction_error": q_pred - q_obs,
        "cu_background_observed": cu_bg,
        "cu_all_observed": cu_all,
        "cu_background_msc_from_q": cu_bg_pred,
        "cu_all_predicted_from_q": cu_all_pred,
        "delta_cu_observed": cu_all - cu_bg,
        "delta_cu_predicted": cu_all_pred - cu_bg_pred,
        "relative_delta_cu_observed": (cu_all - cu_bg) / cu_bg,
        "relative_delta_cu_predicted": (cu_all_pred - cu_bg_pred) / cu_bg_pred,
    }


def load_rows() -> list[dict[str, object]]:
    an_dir = ROOT / "empirical/anopheles_2la/results"
    an_manifest = json.loads((an_dir / "stage4r_manifest.json").read_text())
    an_qc = read_tsv(an_dir / "stage4r_window_qc.tsv")
    an_counts = {k: sum(r["region_class"] == k and r["status"] == "ok" for r in an_qc) for k in ("inside", "outside")}
    assert an_counts == {"inside": 430, "outside": 546}, an_counts
    assert (an_manifest["n_inside_trees"], an_manifest["n_outside_trees"], an_manifest["n_clean_trees"]) == (430, 546, 976)
    an_source = read_tsv(an_dir / "stage4r_fixed_split_scores.tsv")
    an_labels = ["5 taxa | melas + merus", "3 taxa | melas + merus + quadriannulatus", "4 taxa | coluzzii + gambiae"]
    rows: list[dict[str, object]] = []
    for source, display in zip(an_source, an_labels):
        branch_id = source["split"]
        rows.append(row("Anopheles 2La", branch_id, display, 546, 430,
                        float(source["q_baseline_outside"]), float(source["q_baseline_inside"]),
                        float(source["q_baseline_all"]), float(source["CULength_outside"]),
                        float(source["CULength_all"])))

    fire_dir = ROOT / "empirical/fire_ants_chr16/results"
    fire_manifest = json.loads((fire_dir / "stage6c_manifest.json").read_text())
    counts = fire_manifest["raxml_branch_length_audit"]["treatment_counts"]
    assert (counts["T_background"], counts["T_supergene"], counts["T_all"]) == (161, 52, 213), counts
    methods = (fire_dir / "fire_ants_final_methods.md").read_text()
    assert "213 four-BUSCO-gene" in methods and "161-window background" in methods
    cu_source = read_tsv(fire_dir / "stage6c_cu_vs_su_comparison.tsv")
    q_sources = {name: q1_by_split(fire_dir / "astral4_su" / f"T_{name}.fixed_background_su.nwk")
                 for name in ("background", "supergene", "all")}
    fire_labels = {
        "focal_richteri_SB_Sb_pair": "richteri SB/Sb",
        "focal_invicta_macdonaghi_SB_Sb_pair": "invicta/macdonaghi SB/Sb",
        "focal_four_vs_outgroups": "four focal groups | outgroups",
        "focal_related_internal_branch": "geminata | pusillignis",
    }
    for source in cu_source:
        branch_id = parse_split(source["split_id"])
        for name, qmap in q_sources.items():
            if branch_id not in qmap:
                raise ValueError(f"fire-ant split missing from {name}: {branch_id}")
        rows.append(row("Fire ants chr16", branch_id, fire_labels[source["branch_role"]], 161, 52,
                        q_sources["background"][branch_id], q_sources["supergene"][branch_id],
                        q_sources["all"][branch_id], float(source["cu_background"]),
                        float(source["cu_all"])))
    assert len(rows) == 7
    return rows


def write_table(rows: list[dict[str, object]]) -> Path:
    path = RESULTS / "cross_dataset_cu_distortion.tsv"
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: f"{v:.12g}" if isinstance(v, float) else v for k, v in r.items()})
    return path


def write_text(rows: list[dict[str, object]]) -> tuple[Path, Path]:
    caption = RESULTS / "cross_dataset_cu_distortion_caption.md"
    caption.write_text(
        "**Figure. Quartet-mixture prediction of apparent coalescent-unit branch-length distortion across two empirical systems.** "
        "Each point uses a focal branch whose topology is retained when a recombination-suppressed region is added. "
        "For Anopheles 2La, the affected region is the 430-window inside set added to 546 outside windows; for fire ants, "
        "the 52-window chromosome-16 supergene set is added to 161 background windows. The same fixed background split is "
        "scored in all treatments. We calculate pooled support as `q_mix = (1-epsilon) q_bg + epsilon q_R` and transform it as "
        "`t_tilde = -log[3/2 (1-q_mix)]`. (A) Predicted versus observed pooled ASTRAL4 CULength. (B) Predicted versus observed "
        "relative change from the background CULength. Points are colored by dataset; labels identify the biologically focal "
        "branches. The fire-ant data include both strong shortening and a small lengthening example. This is a theory-guided "
        "empirical demonstration of the general quartet-mixture to apparent-CU distortion result, without requiring an exact "
        "2:2 arrangement configuration. It is not independent fitting or validation of all MSRC parameters because `q_R` is "
        "measured from the affected-region data. The fire-ant source study interprets the supergene pattern as adaptive "
        "introgression; this figure does not claim MSRC without gene flow."
    )
    results = RESULTS / "cross_dataset_cu_distortion_results.md"
    max_err = max(abs(float(r["q_prediction_error"])) for r in rows)
    lines = [
        "# Cross-dataset quartet-mixture CU distortion",
        "",
        "Adding a localized genealogy regime changes concordant quartet support while the pooled species-tree topology remains unchanged in both datasets. The linear mixture exactly reconstructs the observed pooled support from the frozen background and affected-region scores, and the MSC quartet-to-CU transform predicts the direction and approximate magnitude of the observed ASTRAL4 CULength change.",
        "",
        "The fire-ant fixed-background analysis provides two strong shortening examples and a small lengthening example. These results support the general branch-length-distortion theorem and do not require an exact 2:2 arrangement configuration. They are a theory-guided empirical demonstration, not independent fitting or validation of all MSRC parameters, because the affected-region support `q_R` is observed rather than predicted from independently estimated MSRC parameters. The source fire-ant study interprets its pattern as adaptive introgression; MSRC without gene flow is not claimed here.",
        "",
        f"Maximum absolute pooled-support prediction error: `{max_err:.3e}`.",
        "",
        "All seven branch transforms were defined (`q_bg` and `q_mix` > 1/3); no branch required an undefined MSC length placeholder.",
        "",
        "| dataset | branch | q_bg | q_R | q_pred | q_obs | CU_bg | CU_all | predicted % change | observed % change |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        label = str(r["display_label"]).replace("|", "/")
        lines.append(f"| {r['dataset']} | {label} | {r['q_background']:.6f} | {r['q_affected']:.6f} | {r['q_all_predicted']:.6f} | {r['q_all_observed']:.6f} | {r['cu_background_observed']:.5f} | {r['cu_all_observed']:.5f} | {100*r['relative_delta_cu_predicted']:.2f}% | {100*r['relative_delta_cu_observed']:.2f}% |")
    results.write_text("\n".join(lines) + "\n")
    return caption, results


def make_figure(rows: list[dict[str, object]]) -> tuple[Path, Path]:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11,
                         "axes.labelsize": 10, "pdf.fonttype": 42, "ps.fonttype": 42})
    colors = {"Anopheles 2La": "#0072B2", "Fire ants chr16": "#D55E00"}
    markers = {"Anopheles 2La": "o", "Fire ants chr16": "s"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), constrained_layout=True)
    ax = axes[0]
    x = [float(r["cu_all_predicted_from_q"]) for r in rows]
    y = [float(r["cu_all_observed"]) for r in rows]
    lim = (0, max(max(x), max(y)) * 1.07)
    ax.plot(lim, lim, color="0.45", lw=1, ls="--", zorder=0)
    for r in rows:
        ax.scatter(r["cu_all_predicted_from_q"], r["cu_all_observed"], s=34, c=colors[r["dataset"]], marker=markers[r["dataset"]], edgecolor="white", linewidth=0.45, zorder=3)
        if r["display_label"] in {"richteri SB/Sb", "invicta/macdonaghi SB/Sb", "5 taxa | melas + merus"}:
            ax.annotate(r["display_label"], (r["cu_all_predicted_from_q"], r["cu_all_observed"]), xytext=(4, 4), textcoords="offset points", fontsize=7, color=colors[r["dataset"]])
    ax.set(xlim=lim, ylim=lim, xlabel="Predicted pooled CU length", ylabel="Observed pooled ASTRAL4 CU length", title="A  Pooled CU length")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.18, lw=0.5)
    ax = axes[1]
    vals = [100 * float(r["relative_delta_cu_predicted"]) for r in rows] + [100 * float(r["relative_delta_cu_observed"]) for r in rows]
    span = max(abs(min(vals)), abs(max(vals))) * 1.14
    ax.plot([-span, span], [-span, span], color="0.45", lw=1, ls="--", zorder=0)
    for r in rows:
        xp = 100 * float(r["relative_delta_cu_predicted"])
        yo = 100 * float(r["relative_delta_cu_observed"])
        ax.scatter(xp, yo, s=34, c=colors[r["dataset"]], marker=markers[r["dataset"]], edgecolor="white", linewidth=0.45, zorder=3)
        if r["display_label"] in {"richteri SB/Sb", "invicta/macdonaghi SB/Sb"}:
            ax.annotate(r["display_label"], (xp, yo), xytext=(4, 4), textcoords="offset points", fontsize=7, color=colors[r["dataset"]])
    ax.axhline(0, color="0.78", lw=0.6); ax.axvline(0, color="0.78", lw=0.6)
    ax.set(xlim=(-span, span), ylim=(-span, span), xlabel="Predicted relative CU change (%)", ylabel="Observed relative CU change (%)", title="B  Relative CU change")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.18, lw=0.5)
    handles = [plt.Line2D([], [], marker=markers[d], color="none", markerfacecolor=colors[d], markeredgecolor="white", markersize=6, label=d) for d in colors]
    axes[1].legend(handles=handles, frameon=False, loc="lower right", handletextpad=0.3, borderpad=0.2)
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"cross_dataset_cu_distortion.{suffix}", dpi=300 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)
    return FIGURES / "cross_dataset_cu_distortion.png", FIGURES / "cross_dataset_cu_distortion.pdf"


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    for r in rows:
        for key in ("q_background", "q_affected", "q_all_predicted", "q_all_observed"):
            assert 0 <= float(r[key]) <= 1, (r["branch_id"], key, r[key])
        assert abs((1 - float(r["epsilon"])) + float(r["epsilon"]) - 1) < 1e-15
        assert abs(float(r["q_prediction_error"])) < 3e-6
    table = write_table(rows)
    caption, results = write_text(rows)
    png, pdf = make_figure(rows)
    print(json.dumps({"table": str(table), "caption": str(caption), "results": str(results), "png": str(png), "pdf": str(pdf), "rows": rows, "max_abs_q_error": max(abs(float(r["q_prediction_error"])) for r in rows)}, indent=2))


if __name__ == "__main__":
    main()
