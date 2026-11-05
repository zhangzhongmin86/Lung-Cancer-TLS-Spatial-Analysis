#!/usr/bin/env python3
"""Intersect NMF communication programs with cell-type co-occurrence support."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


GROUP_MARKERS = ["o", "s", "^", "D", "v", "X", "P", "*"]
GROUP_COLORS = [
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
    "#000000",
]
NODE_COLORS = [
    "#4E79A7",
    "#F28E2B",
    "#59A14F",
    "#E15759",
    "#76B7B2",
    "#B07AA1",
    "#EDC948",
    "#9C755F",
    "#BAB0AC",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modules", required=True, help="nmf_module_interactions.csv")
    parser.add_argument("--cooccurrence-dir", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--top-k", type=int, default=20)
    return parser.parse_args()


def group_from_filename(path: Path) -> str:
    prefix = "specificity_correlation_"
    return path.stem[len(prefix) :].replace("_", " ")


def load_support(directory: Path) -> tuple[dict[frozenset[str], set[str]], list[str]]:
    support: dict[frozenset[str], set[str]] = {}
    groups: list[str] = []
    combined = directory / "cooccurrence_all_groups.csv"
    if combined.exists():
        table = pd.read_csv(combined)
        if "group" in table.columns:
            groups = sorted(table["group"].dropna().astype(str).unique().tolist())
            for row in table.itertuples(index=False):
                pair = frozenset((str(row.Cell_Type_1), str(row.Cell_Type_2)))
                support.setdefault(pair, set()).add(str(row.group))
            return support, groups
    for path in sorted(directory.glob("specificity_correlation_*.csv")):
        group = group_from_filename(path)
        groups.append(group)
        table = pd.read_csv(path)
        if table.empty:
            continue
        for row in table.itertuples(index=False):
            pair = frozenset((str(row.Cell_Type_1), str(row.Cell_Type_2)))
            support.setdefault(pair, set()).add(group)
    return support, groups


def abbreviate(groups: list[str]) -> dict[str, str]:
    preferred = {
        "Normal": "N",
        "Tumor": "T",
        "Autoimmune diseases": "AD",
        "Inflammation": "Inf",
        "Normal adjacent tissue": "NAT",
        "Precancerous condition": "Pre",
        "Tumor metastasis": "Met",
    }
    used: set[str] = set()
    output: dict[str, str] = {}
    for group in groups:
        candidate = preferred.get(group)
        if candidate is None:
            candidate = "".join(word[0] for word in group.split() if word).upper() or "G"
        base = candidate
        suffix = 2
        while candidate in used:
            candidate = f"{base}{suffix}"
            suffix += 1
        output[group] = candidate
        used.add(candidate)
    return output


def celltype_parent_map(metadata_file: Path) -> dict[str, str]:
    metadata = pd.read_csv(metadata_file, index_col=0, low_memory=False)
    required = {"celltype_2", "celltype_3"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
    pairs = metadata.loc[:, ["celltype_3", "celltype_2"]].dropna()
    return (
        pairs.groupby("celltype_3", observed=True)["celltype_2"]
        .agg(lambda values: values.value_counts().index[0])
        .astype(str)
        .to_dict()
    )


def draw_module(
    table: pd.DataFrame,
    module_id: int,
    output_dir: Path,
    parent_map: dict[str, str],
    group_abbr: dict[str, str],
) -> None:
    graph = nx.MultiDiGraph()
    for row in table.itertuples(index=False):
        graph.add_edge(
            str(row.source),
            str(row.target),
            ligand=str(row.ligand),
            receptor=str(row.receptor),
            weight=float(row.weight),
            groups=str(row.supported_groups).split(";"),
        )
    if graph.number_of_edges() == 0:
        return

    parent_types = sorted({parent_map.get(node, "Unknown") for node in graph.nodes})
    parent_colors = {
        parent: NODE_COLORS[index % len(NODE_COLORS)]
        for index, parent in enumerate(parent_types)
    }
    group_colors = {
        group: GROUP_COLORS[index % len(GROUP_COLORS)]
        for index, group in enumerate(group_abbr)
    }
    group_markers = {
        group: GROUP_MARKERS[index % len(GROUP_MARKERS)]
        for index, group in enumerate(group_abbr)
    }

    figure, axis = plt.subplots(figsize=(11, 9))
    positions = nx.circular_layout(graph, scale=0.78)
    node_sizes = []
    node_colors = []
    for node in graph.nodes:
        incident = [data["weight"] for *_, data in graph.edges(node, data=True)]
        incident += [data["weight"] for *_, data in graph.in_edges(node, data=True)]
        node_sizes.append(1500 + 3500 * max(incident, default=0))
        node_colors.append(parent_colors[parent_map.get(node, "Unknown")])

    nx.draw_networkx_nodes(
        graph,
        positions,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors="black",
        linewidths=0.6,
        ax=axis,
    )
    nx.draw_networkx_labels(graph, positions, font_size=9, ax=axis)

    edges = list(graph.edges(keys=True, data=True))
    lr_pairs = sorted({f"{data['ligand']}-{data['receptor']}" for *_, data in edges})
    cmap = plt.get_cmap("tab20")
    lr_colors = {pair: cmap(index % 20) for index, pair in enumerate(lr_pairs)}
    pair_counts: dict[tuple[str, str], int] = {}
    pair_index: dict[tuple[str, str], int] = {}
    for source, target, _, _ in edges:
        pair_counts[(source, target)] = pair_counts.get((source, target), 0) + 1

    for source, target, key, data in edges:
        pair = (source, target)
        index = pair_index.get(pair, 0)
        pair_index[pair] = index + 1
        count = pair_counts[pair]
        curvature = 0.18 * (index - (count - 1) / 2) if count > 1 else 0.05
        lr_pair = f"{data['ligand']}-{data['receptor']}"
        nx.draw_networkx_edges(
            graph,
            positions,
            edgelist=[(source, target, key)],
            edge_color=[lr_colors[lr_pair]],
            width=0.7 + 2.0 * float(data["weight"]),
            arrows=True,
            arrowsize=13,
            arrowstyle="-|>",
            connectionstyle=f"arc3,rad={curvature}",
            min_source_margin=18,
            min_target_margin=18,
            ax=axis,
        )
        start = np.asarray(positions[source], dtype=float)
        end = np.asarray(positions[target], dtype=float)
        midpoint = (start + end) / 2
        direction = end - start
        normal = np.array([-direction[1], direction[0]])
        norm = np.linalg.norm(normal)
        if norm:
            normal /= norm
        groups = [group for group in data["groups"] if group]
        for marker_index, group in enumerate(groups):
            offset = (marker_index - (len(groups) - 1) / 2) * 0.04
            point = midpoint + normal * offset
            axis.scatter(
                point[0],
                point[1],
                s=38,
                marker=group_markers[group],
                c=[group_colors[group]],
                edgecolors="black",
                linewidths=0.35,
                zorder=5,
            )

    parent_legend = [
        mlines.Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markerfacecolor=color,
            markeredgecolor="black",
            label=parent,
        )
        for parent, color in parent_colors.items()
    ]
    group_legend = [
        mlines.Line2D(
            [],
            [],
            marker=group_markers[group],
            linestyle="None",
            markerfacecolor=group_colors[group],
            markeredgecolor="black",
            label=f"{group_abbr[group]}: {group}",
        )
        for group in group_abbr
    ]
    context_legend = axis.legend(
        handles=parent_legend + group_legend,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
        fontsize=8,
    )
    axis.add_artist(context_legend)
    lr_legend = [
        mlines.Line2D([], [], color=lr_colors[pair], linewidth=2, label=pair)
        for pair in lr_pairs
    ]
    axis.legend(
        handles=lr_legend,
        title="Ligand-receptor",
        loc="lower left",
        bbox_to_anchor=(1.01, 0),
        frameon=False,
        fontsize=6,
        title_fontsize=7,
        ncol=2 if len(lr_legend) > 10 else 1,
    )
    axis.set_title(f"COMMUNE module {module_id}")
    axis.set_xlim(-1.02, 1.02)
    axis.set_ylim(-1.02, 1.02)
    axis.axis("off")
    figure.tight_layout()
    for suffix in ("pdf", "png"):
        figure.savefig(
            output_dir / f"module_{module_id}_network.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(figure)


def main() -> None:
    args = parse_args()
    modules_file = Path(args.modules).resolve()
    cooccurrence_dir = Path(args.cooccurrence_dir).resolve()
    metadata_file = Path(args.metadata).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for old_plot in output_dir.glob("module_*_network.*"):
        if old_plot.suffix in {".pdf", ".png"}:
            old_plot.unlink()

    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "sans-serif",
            "axes.unicode_minus": False,
        }
    )
    modules = pd.read_csv(modules_file)
    required = {"module", "weight", "source", "target", "ligand", "receptor"}
    missing = required.difference(modules.columns)
    if missing:
        raise ValueError(f"Module table is missing columns: {sorted(missing)}")

    support, groups = load_support(cooccurrence_dir)
    group_abbr = abbreviate(groups)
    parent_map = celltype_parent_map(metadata_file)
    modules["supported_groups"] = modules.apply(
        lambda row: ";".join(
            sorted(support.get(frozenset((str(row.source), str(row.target))), set()))
        ),
        axis=1,
    )
    constrained = modules.loc[modules["supported_groups"].ne("")].copy()
    constrained = (
        constrained.sort_values(["module", "weight"], ascending=[True, False])
        .groupby("module", observed=True, group_keys=False)
        .head(args.top_k)
    )
    constrained.to_csv(output_dir / "constrained_module_interactions.csv", index=False)

    if constrained.empty:
        summary = pd.DataFrame(columns=["module", "interactions", "cell_types"])
    else:
        summary = constrained.groupby("module", observed=True).apply(
            lambda table: pd.Series(
                {
                    "interactions": len(table),
                    "cell_types": len(set(table["source"]) | set(table["target"])),
                }
            ),
            include_groups=False,
        ).reset_index()
    summary.to_csv(output_dir / "constrained_module_summary.csv", index=False)
    for module_id, table in constrained.groupby("module", observed=True):
        draw_module(table, int(module_id), output_dir, parent_map, group_abbr)

    manifest = {
        "input_nmf_modules": int(modules["module"].nunique()),
        "retained_modules": int(constrained["module"].nunique()),
        "retained_interactions": int(len(constrained)),
        "cooccurrence_groups": groups,
        "top_k_per_module": args.top_k,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
