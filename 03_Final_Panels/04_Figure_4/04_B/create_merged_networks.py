#!/usr/bin/env python3
"""Figure 4B: the five module networks, one above the other, each drawn from
its own Jaccard similarity matrix.

Node names are a lineage letter (defined in the Figure 4 legend) plus the
cluster digit, e.g. 'H9', 'C2', 'P0'; main() checks they are all distinct.
For a byte-identical drawing run with PYTHONHASHSEED=0 (see induced()).

Input:  Module_{1..5}_jaccard_matrix.csv (beside this script)
Output: merged_networks_k5.{svg,pdf,png} in this directory
Run:    PYTHONHASHSEED=0 python create_merged_networks.py
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import networkx as nx
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
import panel_style_cns as style           # noqa: E402
import slots                              # noqa: E402
from shared.node_letters import NODE_LETTER   # noqa: E402

PANEL_LETTER = "B"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

EARLIER_CANVAS_H_MM = 120.0     # edge widths scale by PANEL_H_MM / this

SIMILARITY_THRESHOLD = 0.01
#: Node circle diameter; holds a 6 pt two-character name.
NODE_D_MM = 4.0
#: Axis limit either side of the unit ring; the ring's outer edge is
#: (1 + NODE_D_MM / 2 / unit) and must stay inside it.
LIMIT = 1.35
EDGE_WIDTH_MULTIPLIER = 2
EDGE_ALPHA = 0.7
EDGE_LABEL_THRESHOLD = 1.1  # Effectively disables edge labels (max Jaccard is 1.0)

MODULE_DISPLAY = {1: 'IM-T/NK/DC', 2: 'IM-MoMac', 3: 'IM-Mixed', 4: 'IM-Neutrophil', 5: 'IM-B/Plasma'}
MAX_NODES = 8
DRAWN_LABELS = []       # every node name drawn, for the uniqueness check
LAYOUT_SCALE = 1.0

CELL_TYPE_COLORS = {
    'CD4_T': '#d4a373',
    'CD8_T': '#d4a373',
    'B': '#edafb8',
    'DC': '#dedbd2',
    'MoMac': '#f7e1d7',
    'Neu': '#b0c4b1',
    'NK': '#d4a373',
    'Plasma': '#edafb8',
    'Mast': '#d4a373',
    'Prolif': '#d4a373',
    'Other': '#9E9E9E'
}


def get_cell_type_from_state(cell_state):
    """Extract cell type from cell state name."""
    if 'CD4' in cell_state or 'Th' in cell_state or 'Treg' in cell_state or 'Tfh' in cell_state or 'Tcm' in cell_state or 'Tem' in cell_state or 'Temra' in cell_state:
        return "CD4_T"
    elif 'CD8' in cell_state:
        return "CD8_T"
    elif 'NK' in cell_state:
        return "NK"
    elif 'DC' in cell_state:
        return "DC"
    elif 'Mono' in cell_state or 'Mac' in cell_state or 'MoMac' in cell_state:
        return "MoMac"
    elif 'Neu' in cell_state:
        return "Neu"
    elif cell_state.startswith('C') and '_B_' in cell_state:
        return "B"
    elif 'Plasma' in cell_state:
        return "Plasma"
    elif 'Mast' in cell_state:
        return "Mast"
    elif 'Prolif' in cell_state:
        return "Prolif"
    else:
        return "Other"


def get_short_label(cell_state):
    """One letter for the lineage (shared/node_letters.py) and the cluster
    digit: 'C2_CD8_MAIT_KLRB1' -> 'C2'; 'C9_CD4_Th17_IL17A' -> 'H9'.
    """
    parts = cell_state.split('_')
    digit = parts[0][1:]
    lineage = parts[1]
    if lineage not in NODE_LETTER:
        raise KeyError(f"{cell_state}: no letter for lineage {lineage!r} in "
                       f"labels.NODE_LETTER")
    return f"{NODE_LETTER[lineage]}{digit}"


def induced(G, nodes):
    """The subgraph on `nodes`, in the order `G` itself holds them.

    `Graph.subgraph` filters through a set, so node order (and the circular
    layout) would vary between runs. Edge paint order still depends on
    PYTHONHASHSEED.
    """
    keep = set(nodes)
    sub = G.subgraph(keep)
    H = nx.Graph()
    H.add_nodes_from((n, dict(sub.nodes[n])) for n in G.nodes() if n in keep)
    H.add_edges_from(sub.edges(data=True))
    return H


def create_network_graph(jaccard_df, similarity_threshold):
    """Create NetworkX graph from Jaccard similarity matrix."""
    G = nx.Graph()

    for state in jaccard_df.index:
        G.add_node(state)

    for i, state1 in enumerate(jaccard_df.index):
        for j, state2 in enumerate(jaccard_df.columns):
            if i < j:  # Avoid duplicates
                similarity = jaccard_df.loc[state1, state2]
                if similarity >= similarity_threshold:
                    G.add_edge(state1, state2, weight=similarity)

    return G


def plot_single_network(ax, G, module_id, mark):
    """Plot a single network on a given axis."""
    if len(G.nodes()) > MAX_NODES:
        degrees = dict(G.degree())
        top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:MAX_NODES]
        G = induced(G, top_nodes)

    if len(G.nodes()) == 0:
        ax.text(0.5, 0.5, f'{MODULE_DISPLAY[module_id]}\n(No nodes)',
                ha='center', va='center', fontsize=style.body_pt())
        ax.axis('off')
        return

    pos = nx.circular_layout(G)
    # Scale positions to bring nodes closer together
    pos = {node: (x * LAYOUT_SCALE, y * LAYOUT_SCALE) for node, (x, y) in pos.items()}

    node_colors = [CELL_TYPE_COLORS.get(get_cell_type_from_state(node), '#CCCCCC')
                   for node in G.nodes()]

    nodes = nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=(NODE_D_MM * style.PT_PER_MM) ** 2,   # s = diameter^2 pt
        alpha=0.9,
        edgecolors='black',
        linewidths=style.EDGE_PT,
        ax=ax
    )
    nodes.set_clip_on(False)

    edges = G.edges()
    if len(edges) > 0:
        weights = [G[u][v]['weight'] for u, v in edges]
        edge_widths = [w * EDGE_WIDTH_MULTIPLIER * mark for w in weights]

        nx.draw_networkx_edges(
            G, pos,
            width=edge_widths,
            alpha=EDGE_ALPHA,
            edge_color='gray',
            ax=ax
        )

        edge_labels = {}
        for u, v in edges:
            weight = G[u][v]['weight']
            if weight >= EDGE_LABEL_THRESHOLD:
                edge_labels[(u, v)] = f"{weight:.2f}"

        if edge_labels:
            nx.draw_networkx_edge_labels(
                G, pos,
                edge_labels=edge_labels,
                font_size=style.tick_pt(),
                ax=ax
            )

    labels = {node: get_short_label(node) for node in G.nodes()}
    nx.draw_networkx_labels(
        G, pos,
        labels=labels,
        font_size=style.tick_pt(),
        ax=ax
    )
    DRAWN_LABELS.extend(labels.values())

    ax.set_title(MODULE_DISPLAY[module_id], pad=3 * mark)

    ax.set_xlim(-LIMIT, LIMIT)
    ax.set_ylim(-LIMIT, LIMIT)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')


def main():
    """Generate merged network figure."""
    script_dir = Path(__file__).parent

    style.apply(title_fontsize=7, fontsize_legend=6, legend_fontsize=6)
    MARK = PANEL_H_MM / EARLIER_CANVAS_H_MM

    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 5, 1)

    for module_id in range(1, 6):
        jaccard_path = script_dir / f'Module_{module_id}_jaccard_matrix.csv'
        jaccard_df = pd.read_csv(jaccard_path, index_col=0)
        G = create_network_graph(jaccard_df, SIMILARITY_THRESHOLD)
        plot_single_network(axes[module_id - 1], G, module_id, MARK)

    if len(set(DRAWN_LABELS)) != len(DRAWN_LABELS):
        dup = sorted({n for n in DRAWN_LABELS if DRAWN_LABELS.count(n) > 1})
        raise RuntimeError(f"node names are not unique: {dup}")

    style.fit_margins(fig, pad_mm=0.6, cell_mm=LETTER_CELL)
    over = style.overflow_mm(fig)
    if max(over) > 0:
        raise RuntimeError(
            f"ink outside the {PANEL_W_MM} x {PANEL_H_MM} mm canvas "
            f"(l,r,b,t mm): {over}")
    intruders = style.letter_clear(fig, LETTER_CELL)
    if intruders:
        raise RuntimeError(f"ink under the panel letter cell: {intruders}")

    stem = 'merged_networks_k5'
    style.save_panel(fig, script_dir / stem)
    print(f"  Saved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm; "
          f"{len(DRAWN_LABELS)} node names, all distinct")


if __name__ == "__main__":
    main()
