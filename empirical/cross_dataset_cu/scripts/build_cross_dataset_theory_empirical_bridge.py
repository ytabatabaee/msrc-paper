#!/usr/bin/env python3
"""Build the frozen three-system theory-to-empirical bridge figure."""

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
    a, b = tuple(sorted(left)), tuple(sorted(right))
    chosen, other = (a, b) if (len(a), a) <= (len(b), b) else (b, a)
    return ",".join(chosen) + "|" + ",".join(other)


def parse_split(text: str) -> str:
    left, right = text.split("|", 1)
    return canonical_split(set(left.split(",")), set(left.split(",")) | set(right.split(",")))


class NewickParser:
    def __init__(self, text: str):
        self.text, self.i = text.strip().rstrip(";"), 0
        self.records: list[tuple[set[str], dict[str, float]]] = []

    def parse(self) -> set[str]:
        leaves = self.subtree()
        self.skip_space()
        if self.i != len(self.text):
            raise ValueError(f"unparsed Newick suffix at {self.i}")
        return leaves

    def skip_space(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def skip_length(self) -> None:
        self.skip_space()
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            while self.i < len(self.text) and self.text[self.i] not in ",()":
                self.i += 1

    def subtree(self) -> set[str]:
        self.skip_space()
        if self.text[self.i] == "(":
            self.i += 1
            leaves: set[str] = set()
            while True:
                leaves |= self.subtree()
                self.skip_space()
                if self.text[self.i] == ",":
                    self.i += 1
                elif self.text[self.i] == ")":
                    self.i += 1
                    break
                else:
                    raise ValueError(f"expected comma or close at {self.i}")
            self.skip_space()
            annotation: dict[str, float] = {}
            if self.i < len(self.text) and self.text[self.i] == "'":
                self.i += 1
                start = self.i
                end = self.text.index("'", start)
                annotation = {key: float(value.strip("[]")) for key, value in re.findall(
                    r"([A-Za-z][A-Za-z0-9_]*)=([^;]+)", self.text[start:end]
                )}
                self.i = end + 1
            self.skip_length()
            if annotation:
                self.records.append((leaves, annotation))
            return leaves
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",():;'":
            self.i += 1
        label = self.text[start:self.i].strip()
        if not label:
            raise ValueError(f"empty leaf at {self.i}")
        self.skip_length()
        return {label}


def q1_by_split(path: Path) -> dict[str, tuple[float, float, float]]:
    parser = NewickParser(path.read_text())
    all_taxa = parser.parse()
    values: dict[str, tuple[float, float, float]] = {}
    for leaves, annotation in parser.records:
        if "q1" not in annotation:
            continue
        key = canonical_split(leaves, all_taxa)
        if key in values:
            raise ValueError(f"duplicate split in {path}: {key}")
        values[key] = (annotation["q1"], annotation["q2"], annotation["q3"])
    return values


def threshold(b_s: float, b_j: float, r_s: float, r_j: float) -> float | None:
    if not (b_s > b_j and r_j > r_s):
        return None
    value = (b_s - b_j) / ((b_s - b_j) + (r_j - r_s))
    if not 0 <= value <= 1:
        raise AssertionError(value)
    return value


def mixture(b: tuple[float, float, float], r: tuple[float, float, float], epsilon: float) -> tuple[float, float, float]:
    return tuple((1 - epsilon) * x + epsilon * y for x, y in zip(b, r))


def topology(q: tuple[float, float, float], names=("Q_SPECIES", "Q_ALT", "Q_OTHER")) -> str:
    return names[max(range(3), key=lambda i: q[i])]


def make_threshold_row(dataset: str, branch_id: str, display: str, definition: str, epsilon: float,
                       b: tuple[float, float, float], r: tuple[float, float, float],
                       pooled: tuple[float, float, float], background_topology: str,
                       pooled_topology: str, affected_topology: str, notes: str) -> dict[str, object]:
    stars = [threshold(b[0], b[i], r[0], r[i]) for i in (1, 2)]
    valid = [(i, x) for i, x in zip((1, 2), stars) if x is not None]
    if not valid:
        raise ValueError(f"finite threshold required for included row: {branch_id}")
    winner, epsilon_star = min(valid, key=lambda item: item[1])
    return {
        "dataset": dataset,
        "branch_id": branch_id,
        "display_label": display,
        "epsilon_definition": definition,
        "epsilon_observed": epsilon,
        "qS_background": b[0],
        "qA_background": b[1],
        "qO_background": b[2],
        "qS_affected": r[0],
        "qA_affected": r[1],
        "qO_affected": r[2],
        "winning_affected_alternative": "Q_ALT" if winner == 1 else "Q_OTHER",
        "epsilon_star_alt1": stars[0],
        "epsilon_star_alt2": stars[1],
        "epsilon_star": epsilon_star,
        "observed_below_or_above_threshold": "below" if epsilon < epsilon_star else "above",
        "background_topology": background_topology,
        "pooled_topology": pooled_topology,
        "affected_only_topology": affected_topology,
        "notes": notes,
    }


def load_mouse() -> tuple[dict[str, object], dict[str, float]]:
    result_dir = ROOT / "empirical/house_mouse_t_complex/results"
    contrib = read_tsv(result_dir / "stage2_arrangement_pattern_contributions.tsv")
    patterns = {r["status_pattern"]: r for r in contrib}
    assert set(patterns) == {"SSS", "SST", "STS", "STT", "TSS", "TST", "TTS", "TTT"}
    homogeneous = ("SSS", "TTT")
    mixed = ("SST", "STS", "STT", "TSS", "TST", "TTS")

    def normalized(names: tuple[str, ...]) -> tuple[float, float, float]:
        weights = [float(patterns[n]["pattern_weight"]) for n in names]
        total = sum(weights)
        vals = [sum(float(patterns[n][f"contribution_q_{q}"]) for n in names) / total for q in ("species", "t_alt", "other")]
        # The frozen quartet summaries retain a small unresolved component;
        # q_B and q_R therefore sum to the resolved fraction, as in the source.
        assert 0 < sum(vals) <= 1 + 1e-10
        return tuple(vals)

    q_b = normalized(homogeneous)
    q_r = normalized(mixed)
    epsilon = sum(float(patterns[n]["pattern_weight"]) for n in mixed)
    assert abs(epsilon - (1 - (0.325 + 0.02857142857))) < 1e-8
    all_summary = next(r for r in read_tsv(result_dir / "stage2_fixed_quartet_summary.tsv") if r["treatment"] == "ALL_TIPS")
    q_all = tuple(float(all_summary[f"mean_q_{q}"]) for q in ("species", "t_alt", "other"))
    reconstructed = mixture(q_b, q_r, epsilon)
    error = max(abs(x - y) for x, y in zip(reconstructed, q_all))
    assert error < 2e-8, (reconstructed, q_all, error)
    astral = read_tsv(result_dir / "stage2_astral_summary.tsv")
    topo = {r["treatment"]: r["topology"] for r in astral}
    assert topo["T0_STANDARD"] == "Q_SPECIES" and topo["T1_ALL_WINDOWS"] == "Q_T_ALT"
    row = make_threshold_row(
        "House mouse chr17 t-complex", "mixed_arrangement_quartets", "mixed S/T quartets",
        "weight of SST, STS, STT, TSS, TST, and TTS induced quartets", epsilon, q_b, q_r,
        q_all, "Q_SPECIES", "Q_T_ALT", topology(q_r),
        "SSS and TTT are homogeneous; mixed classes drive the switch.",
    )
    row["winning_affected_alternative"] = "Q_T_ALT"
    row["affected_only_topology"] = "Q_T_ALT"
    return {"row": row, "q_b": q_b, "q_r": q_r, "q_all": q_all, "reconstructed": reconstructed, "max_error": error}, patterns


def load_threshold_rows() -> tuple[list[dict[str, object]], dict[str, object], list[str]]:
    rows: list[dict[str, object]] = []
    excluded: list[str] = []
    an_dir = ROOT / "empirical/anopheles_2la/results"
    manifest = json.loads((an_dir / "stage4r_manifest.json").read_text())
    qc = read_tsv(an_dir / "stage4r_window_qc.tsv")
    counts = {k: sum(r["region_class"] == k and r["status"] == "ok" for r in qc) for k in ("inside", "outside")}
    assert counts == {"inside": 430, "outside": 546}
    assert (manifest["n_inside_trees"], manifest["n_outside_trees"], manifest["n_clean_trees"]) == (430, 546, 976)
    an_topology = {"background": "Q_SPECIES", "pooled": "Q_SPECIES", "affected": "Q_ALT"}
    for source in read_tsv(an_dir / "stage4r_fixed_split_scores.tsv"):
        b = tuple(float(source[f"q_{x}_outside"]) for x in ("baseline", "alt1", "alt2"))
        r = tuple(float(source[f"q_{x}_inside"]) for x in ("baseline", "alt1", "alt2"))
        pooled = tuple(float(source[f"q_{x}_all"]) for x in ("baseline", "alt1", "alt2"))
        stars = [threshold(b[0], b[i], r[0], r[i]) for i in (1, 2)]
        if all(x is None for x in stars):
            excluded.append(f"Anopheles 2La: {source['split']} — neither affected alternative exceeds r_S")
            continue
        affected_top = topology(r)
        rows.append(make_threshold_row(
            "Anopheles 2La", source["split"],
            {"arabiensis,coluzzii,gambiae,quadriannulatus|melas,merus": "5 taxa | melas + merus",
             "arabiensis,melas,merus,quadriannulatus|coluzzii,gambiae": "4 taxa | coluzzii + gambiae"}[source["split"]],
            "inside 2La windows / all usable non-boundary windows", 430 / 976, b, r, pooled,
            an_topology["background"], an_topology["pooled"], affected_top,
            "T_outside6 vs T_all6 RF = 0; T_inside6 differs from T_outside6 (RF = 4).",
        ))

    fire_dir = ROOT / "empirical/fire_ants_chr16/results"
    manifest = json.loads((fire_dir / "stage6c_manifest.json").read_text())
    counts = manifest["raxml_branch_length_audit"]["treatment_counts"]
    assert (counts["T_background"], counts["T_supergene"], counts["T_all"]) == (161, 52, 213)
    free = {r["treatment"]: r for r in read_tsv(fire_dir / "stage6c_free_topology_summary.tsv")}
    assert free["T_background"]["unrooted_RF_to_free_background"] == "0"
    assert free["T_all"]["unrooted_RF_to_free_background"] == "0"
    assert free["T_supergene"]["unrooted_RF_to_free_background"] == "4"
    qmaps = {name: q1_by_split(fire_dir / "astral4_su" / f"T_{name}.fixed_background_su.nwk") for name in ("background", "supergene", "all")}
    labels = {"focal_richteri_SB_Sb_pair": "richteri SB/Sb", "focal_invicta_macdonaghi_SB_Sb_pair": "invicta/macdonaghi SB/Sb"}
    for source in read_tsv(fire_dir / "stage6c_cu_vs_su_comparison.tsv"):
        if source["branch_role"] not in labels:
            excluded.append(f"Fire ants chr16: {source['branch_role']} — not one of the topology-sensitive focal branches requested for Panel B")
            continue
        key = parse_split(source["split_id"])
        b, r, pooled = qmaps["background"][key], qmaps["supergene"][key], qmaps["all"][key]
        rows.append(make_threshold_row(
            "Fire ants chr16", key, labels[source["branch_role"]],
            "supergene windows / all chromosome-16 windows", 52 / 213, b, r, pooled,
            "Q_SPECIES", "Q_SPECIES", topology(r),
            "T_background vs T_all RF = 0; T_supergene vs background RF = 4.",
        ))

    mouse, _ = load_mouse()
    rows.append(mouse["row"])
    assert len(rows) == 5
    return rows, mouse, excluded


def write_threshold_table(rows: list[dict[str, object]]) -> Path:
    path = RESULTS / "cross_dataset_topology_threshold.tsv"
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("NA" if v is None else f"{v:.12g}" if isinstance(v, float) else v) for k, v in row.items()})
    return path


def build_figure(rows: list[dict[str, object]], mouse: dict[str, object]) -> tuple[Path, Path]:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 10,
                         "axes.labelsize": 9, "pdf.fonttype": 42, "ps.fonttype": 42})
    colors = {"Anopheles 2La": "#0072B2", "Fire ants chr16": "#D55E00", "House mouse chr17 t-complex": "#009E73"}
    markers = {"Anopheles 2La": "o", "Fire ants chr16": "s", "House mouse chr17 t-complex": "D"}
    cu_rows = read_tsv(RESULTS / "cross_dataset_cu_distortion.tsv")
    fig = plt.figure(figsize=(8.0, 7.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1, 1.12), hspace=0.35, wspace=0.28)

    ax = fig.add_subplot(grid[0, 0])
    values = [100 * float(r["relative_delta_cu_predicted"]) for r in cu_rows] + [100 * float(r["relative_delta_cu_observed"]) for r in cu_rows]
    lim = max(abs(min(values)), abs(max(values))) * 1.12
    ax.plot([-lim, lim], [-lim, lim], "--", color="0.45", lw=0.9)
    for r in cu_rows:
        x, y = 100 * float(r["relative_delta_cu_predicted"]), 100 * float(r["relative_delta_cu_observed"])
        ax.scatter(x, y, s=28, c=colors[r["dataset"]], marker=markers[r["dataset"]], edgecolor="white", linewidth=0.4, zorder=3)
        if r["display_label"] in {"5 taxa | melas + merus", "richteri SB/Sb", "invicta/macdonaghi SB/Sb"}:
            offsets = {"5 taxa | melas + merus": (4, -11), "richteri SB/Sb": (4, 5), "invicta/macdonaghi SB/Sb": (4, 5)}
            ax.annotate(r["display_label"].replace(" | ", " / "), (x, y), xytext=offsets[r["display_label"]], textcoords="offset points", fontsize=6.8, color=colors[r["dataset"]])
    ax.axhline(0, color="0.82", lw=0.55); ax.axvline(0, color="0.82", lw=0.55)
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="Predicted relative CU change (%)", ylabel="Observed relative CU change (%)", title="A. Quartet mixture predicts CU distortion")
    ax.set_aspect("equal", adjustable="box"); ax.grid(alpha=0.18, lw=0.45)

    ax = fig.add_subplot(grid[0, 1])
    ax.axvspan(0, 0.5, color="#56B4E9", alpha=0.08, zorder=0)
    ax.axvspan(0.5, 1, color="#E69F00", alpha=0.07, zorder=0)
    y_positions = list(range(len(rows)))[::-1]
    for y, r in zip(y_positions, rows):
        e_obs, e_star = float(r["epsilon_observed"]), float(r["epsilon_star"])
        color = colors[r["dataset"]]
        ax.plot([e_obs, e_star], [y, y], color="0.60", lw=1.1, zorder=1)
        ax.scatter(e_obs, y, s=34, c=color, marker="o", edgecolor="white", linewidth=0.4, zorder=3)
        ax.scatter(e_star, y, s=37, c="white", marker="D", edgecolor=color, linewidth=1.2, zorder=4)
        ax.text(1.01, y, "below" if e_obs < e_star else "above", transform=ax.get_yaxis_transform(), va="center", fontsize=6.8, color="#0072B2" if e_obs < e_star else "#D55E00")
    ax.set(yticks=y_positions, yticklabels=[str(r["display_label"]) for r in rows], xlim=(0, 1), ylim=(-1, len(rows)), xlabel="Affected fraction, ε", title="B. Observed fraction versus topology threshold ε*")
    ax.grid(axis="x", alpha=0.2, lw=0.45)
    ax.text(0.25, 0.98, "topology preserved", transform=ax.transAxes, ha="center", va="top", fontsize=7, color="#0072B2")
    ax.text(0.75, 0.98, "topology failure", transform=ax.transAxes, ha="center", va="top", fontsize=7, color="#D55E00")
    ax.legend([plt.Line2D([], [], marker="o", color="none", markerfacecolor="0.55", markersize=5), plt.Line2D([], [], marker="D", color="#444", markerfacecolor="white", markersize=5)], ["ε_obs", "ε*"], frameon=False, loc="lower right", fontsize=7, ncol=2, handletextpad=0.2, columnspacing=0.6)

    ax = fig.add_subplot(grid[1, :])
    q_b, q_r = mouse["q_b"], mouse["q_r"]
    epsilon_star = float(mouse["row"]["epsilon_star"]); epsilon_obs = float(mouse["row"]["epsilon_observed"])
    xs = [i / 500 for i in range(501)]
    lines = {"Q_SPECIES": ("#0072B2", 0), "Q_T_ALT": ("#D55E00", 1), "Q_OTHER": ("#666666", 2)}
    for name, (color, i) in lines.items():
        ys = [(1 - x) * q_b[i] + x * q_r[i] for x in xs]
        ax.plot(xs, ys, color=color, lw=2, label=name)
    q_star = mixture(q_b, q_r, epsilon_star)
    q_obs = tuple(mouse["q_all"])
    ax.axvline(epsilon_star, color="#444", ls="--", lw=0.9)
    ax.axvline(epsilon_obs, color="#009E73", ls=":", lw=1.2)
    ax.scatter([epsilon_obs] * 3, q_obs, c=[lines[n][0] for n in lines], s=32, edgecolor="white", linewidth=0.5, zorder=4)
    ax.text(epsilon_star, 0.505, "ε* = %.3f" % epsilon_star, ha="center", va="bottom", fontsize=7, color="#444")
    ax.text(epsilon_obs, 0.505, "ε_obs = %.3f" % epsilon_obs, ha="center", va="bottom", fontsize=7, color="#009E73")
    ax.text(0.15, 0.47, "Q_SPECIES dominant", color="#0072B2", fontsize=8)
    ax.text(0.72, 0.47, "Q_T_ALT dominant", color="#D55E00", fontsize=8)
    ax.set(xlim=(0, 1), ylim=(0.25, 0.52), xlabel="Mixed-arrangement weight ε", ylabel="Expected quartet frequency", title="C. Mixed arrangements cross the topology threshold in house mouse")
    ax.grid(alpha=0.18, lw=0.45); ax.legend(frameon=False, loc="lower right", ncol=3, fontsize=7)
    ax.annotate("ALL_TIPS", (epsilon_obs, q_obs[1]), xytext=(5, 4), textcoords="offset points", fontsize=7, color="#009E73")

    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"cross_dataset_theory_empirical_bridge.{suffix}", dpi=300 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)
    return FIGURES / "cross_dataset_theory_empirical_bridge.png", FIGURES / "cross_dataset_theory_empirical_bridge.pdf"


def write_text(rows: list[dict[str, object]], mouse: dict[str, object], excluded: list[str]) -> tuple[Path, Path]:
    caption = RESULTS / "cross_dataset_theory_empirical_bridge_caption.md"
    caption.write_text(
        "**Figure. A theory-to-empirical bridge from quartet mixtures to CU distortion and topology failure.** "
        "(A) Predicted versus observed relative CU change for the frozen Anopheles 2La and fire-ant chromosome-16 branches. "
        "The identity-line agreement illustrates the below-threshold regime: pooled topology is retained while apparent CU lengths change. "
        "(B) Observed affected fraction (filled circles) versus the first valid alternative-topology crossing threshold ε* (open diamonds). "
        "Anopheles and fire ants lie below their thresholds, whereas the house-mouse mixed-arrangement fraction lies above it. "
        "(C) For house mouse, the homogeneous SSS/TTT quartet distribution is mixed with the six mixed S/T arrangement classes. "
        "The species topology dominates below ε*, while Q_T_ALT dominates above ε*; the observed ALL_TIPS mixture is marked. "
        "The mouse switch is associated with mixed arrangement classes: homogeneous SSS and TTT each favor Q_SPECIES. "
        "These are empirical genealogy regimes consistent with the general mixture theorem, not claims that the affected component is generated solely by MSRC. "
        "Affected-component distributions are measured from the data, so the figure is a theory-guided empirical demonstration rather than full parameter-level MSRC validation. "
        "The fire-ant source study interprets its pattern as adaptive introgression, and Anopheles has known introgression; neither system is presented as proof of a rearrangement-only mechanism."
    )
    results = RESULTS / "cross_dataset_theory_empirical_bridge_results.md"
    lines = [
        "# Cross-dataset theory–empirical bridge",
        "",
        "The quartet mixture predicts two linked consequences. When the affected fraction is below the first valid crossing threshold ε*, the pooled topology can remain stable while the concordant support and apparent CU branch length change. When ε exceeds ε*, an alternative topology overtakes the background topology.",
        "",
        "Anopheles 2La uses 430 inside windows among 976 usable non-boundary windows (ε = 0.440574). Both included topology-sensitive branches are below threshold, consistent with `T_outside6` versus `T_all6` RF = 0, while the inside-only topology differs from the outside topology (RF = 4). Fire ants use 52 supergene windows among 213 windows (ε = 0.244131). Both focal branches are below threshold; free-topology Stage 6C gives RF = 0 for background versus all and RF = 4 for supergene versus background.",
        "",
        "House mouse defines ε as the weight of mixed arrangement-state quartets SST, STS, STT, TSS, TST, and TTS within the t-complex, rather than as an inversion-wide fraction. The observed mixed weight is ε = 0.646429, above ε* = 0.311776. The mixture reconstructs the frozen ALL_TIPS vector and agrees with the observed ASTRAL switch from Q_SPECIES (`T0_STANDARD`) to Q_T_ALT (`T1_ALL_WINDOWS`). SSS and TTT are both Q_SPECIES-favoring; the switch is driven by the composition of mixed classes.",
        "",
        "This figure is a theory-guided empirical demonstration, not full parameter-level MSRC validation: the affected-component quartet distributions are measured from the data rather than predicted from independently estimated MSRC parameters. The fire-ant source study’s adaptive-introgression interpretation is retained, and Anopheles has known introgression; neither result establishes a rearrangement-only mechanism.",
        "",
        "| dataset | branch/configuration | ε_obs | ε*_alt1 | ε*_alt2 | ε* | relation | background | pooled | affected-only |",
        "|---|---|---:|---:|---:|---:|---|---|---|---|",
    ]
    for r in rows:
        fmt = lambda x: "NA" if x is None else f"{float(x):.6f}"
        lines.append(f"| {r['dataset']} | {r['display_label']} | {float(r['epsilon_observed']):.6f} | {fmt(r['epsilon_star_alt1'])} | {fmt(r['epsilon_star_alt2'])} | {float(r['epsilon_star']):.6f} | {r['observed_below_or_above_threshold']} | {r['background_topology']} | {r['pooled_topology']} | {r['affected_only_topology']} |")
    lines += ["", f"Mouse reconstructed ALL_TIPS q = ({mouse['reconstructed'][0]:.10f}, {mouse['reconstructed'][1]:.10f}, {mouse['reconstructed'][2]:.10f}); maximum absolute error = {mouse['max_error']:.3e}.", "", "Excluded branches:"]
    lines.extend(f"- {item}" for item in excluded)
    results.write_text("\n".join(lines) + "\n")
    return caption, results


def validate(rows: list[dict[str, object]], mouse: dict[str, object]) -> None:
    assert abs(threshold(0.8, 0.2, 0.3, 0.5) - 0.75) < 1e-12
    assert threshold(0.8, 0.2, 0.5, 0.3) is None
    assert abs(sum(float(mouse["patterns"][x]["pattern_weight"]) for x in ("SSS", "TTT")) - 0.35357142857) < 1e-8
    assert abs(sum(float(mouse["patterns"][x]["pattern_weight"]) for x in ("SST", "STS", "STT", "TSS", "TST", "TTS")) - float(mouse["row"]["epsilon_observed"])) < 1e-8
    assert all(0 <= float(r["epsilon_star"]) <= 1 for r in rows)
    for r in rows:
        assert (r["observed_below_or_above_threshold"] == "below") == (float(r["epsilon_observed"]) < float(r["epsilon_star"]))
    assert sum(r["observed_below_or_above_threshold"] == "below" for r in rows if r["dataset"] != "House mouse chr17 t-complex") == 4
    assert rows[-1]["observed_below_or_above_threshold"] == "above"


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    rows, mouse, excluded = load_threshold_rows()
    mouse["patterns"] = load_mouse()[1]
    validate(rows, mouse)
    table = write_threshold_table(rows)
    png, pdf = build_figure(rows, mouse)
    caption, results = write_text(rows, mouse, excluded)
    print(json.dumps({"table": str(table), "png": str(png), "pdf": str(pdf), "caption": str(caption), "results": str(results), "rows": rows, "mouse_reconstructed_all_tips": mouse["reconstructed"], "mouse_max_error": mouse["max_error"], "excluded": excluded}, indent=2))


if __name__ == "__main__":
    main()
