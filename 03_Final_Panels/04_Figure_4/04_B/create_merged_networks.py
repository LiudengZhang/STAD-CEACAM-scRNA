#!/usr/bin/env python3
"""
Figure 4 panel B - the five module networks, one above the other, each drawn
from its own Jaccard similarity matrix.

  printed panel  Figure 4 B       (PROVENANCE.csv - the directory is "04_B";
                                   do NOT read the directory as the letter)

The five matrices, the similarity threshold, the node selection, the circular
layout, the cell-type colouring, the edge weights and the short node names are
unchanged. Only the canvas and the type change: the panel is drawn at the
millimetre rectangle it prints in and set in the figure's one type system.

THE NODE NAMES ARE AN EXEMPTION
    Eight node names are set around a circle in a row 18.9 mm tall, five rows
    in a column 21.8 mm wide. At the figure's 6 pt floor they collide: five
    overlapping pairs, the worst 4.92 mm2, and five names with glyphs running
    off the drawing. The footprint is the published one and cannot grow -
    measured by redrawing at a range of canvas sizes, the collisions clear only
    at 1.55x it, 33.8 x 146.2 mm - and no abbreviation helps, because dropping
    the cluster prefix leaves four nodes of one module reading Mac and two
    reading Mono. So the size is measured rather than chosen: the panel was
    redrawn from 6.00 pt down in 0.25 pt steps and then refined in 0.05 pt
    steps, and NODE_LABEL_PT is the largest size at which no two names overlap
    and no glyph is clipped or painted over. It is 2.1x the size the published
    page sets these names at. The five module titles are at the figure's own
    body size and are not part of the exemption.

MARK
    The earlier drawing set no SCALE. It drew a 70 x 120 mm canvas of five
    stacked axes, each of which `set_aspect('equal')` reduced to the height of
    its own row - 24 mm, the smaller of the two. The panel prints 94.3 mm tall,
    so a row is 18.9 mm and

        MARK = PANEL_H_MM / EARLIER_CANVAS_H_MM
        AREA = MARK ** 2

    The node marker is an area and is scaled by AREA; the node border, the edge
    widths and the title pad are lengths in points and are scaled by MARK. The
    layout scale and the axis limits are dimensionless and are unchanged, so
    every node sits where it sat.

Input : Module_{1..5}_jaccard_matrix.csv (beside this script)
Output: this directory / merged_networks_k5.{svg,pdf,png}
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
import panel_style_cns as style           # noqa: E402
import slots                              # noqa: E402

PANEL_LETTER = "B"
PANEL_W_MM, PANEL_H_MM = slots.size_mm(4, PANEL_LETTER)
LETTER_CELL = slots.letter_cell_mm(4, PANEL_LETTER)

EARLIER_CANVAS_H_MM = 120.0     # the earlier canvas height, for MARK only

#: The size the eight node names are set at. This is Figure 4's second author
#: exemption and the number is a measurement, not a choice: see THE NODE NAMES
#: ARE AN EXEMPTION in the module docstring.
NODE_LABEL_PT = 3.55

# Parameters
SIMILARITY_THRESHOLD = 0.01
NODE_SIZE = 150  # Smaller nodes
EDGE_WIDTH_MULTIPLIER = 2
EDGE_ALPHA = 0.7
EDGE_LABEL_THRESHOLD = 1.1  # Effectively disables edge labels (max Jaccard is 1.0)

# Module display names
MODULE_DISPLAY = {1: 'IM-T/NK/DC', 2: 'IM-MoMac', 3: 'IM-Mixed', 4: 'IM-Neutrophil', 5: 'IM-B/Plasma'}
MAX_NODES = 8
INCLUDE_STATS = True
LAYOUT_SCALE = 1.0

# Cell type colors (reference style)
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
    """Get shortened label for display."""
    parts = cell_state.split('_')
    if len(parts) >= 3:
        short = f"{parts[0]}_{parts[1]}"
    else:
        short = cell_state

    # Custom label replacements
    label_mappings = {
        'C2_MoMac': 'C2_Mac',
        'C0_Plasma': 'C0_PB'
    }

    return label_mappings.get(short, short)


def induced(G, nodes):
    """The subgraph on `nodes`, in the order `G` itself holds them.

    `Graph.subgraph` keeps its node filter in a set, and a set of strings
    orders by a hash Python randomises per process, so the copy comes back in a
    different order on every run. `circular_layout` seats nodes in that order,
    so the same eight names took a different seat on the circle each time the
    panel was drawn - and a figure that redraws differently on every run cannot
    be reproduced from its deposited code. The nodes, the edges and the weights
    are identical either way; only the seating moved. It is pinned here to the
    order the module's similarity matrix lists them in.
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

    # Add nodes
    for state in jaccard_df.index:
        G.add_node(state)

    # Add edges for similarities above threshold
    for i, state1 in enumerate(jaccard_df.index):
        for j, state2 in enumerate(jaccard_df.columns):
            if i < j:  # Avoid duplicates
                similarity = jaccard_df.loc[state1, state2]
                if similarity >= similarity_threshold:
                    G.add_edge(state1, state2, weight=similarity)

    return G


def plot_single_network(ax, G, module_id, mark, area):
    """Plot a single network on a given axis."""
    # Limit nodes if needed
    if len(G.nodes()) > MAX_NODES:
        degrees = dict(G.degree())
        top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:MAX_NODES]
        G = induced(G, top_nodes)

    if len(G.nodes()) == 0:
        ax.text(0.5, 0.5, f'{MODULE_DISPLAY[module_id]}\n(No nodes)',
                ha='center', va='center', fontsize=style.body_pt())
        ax.axis('off')
        return

    # Circular layout
    pos = nx.circular_layout(G)
    # Scale positions to bring nodes closer together
    pos = {node: (x * LAYOUT_SCALE, y * LAYOUT_SCALE) for node, (x, y) in pos.items()}

    # Get node colors
    node_colors = [CELL_TYPE_COLORS.get(get_cell_type_from_state(node), '#CCCCCC')
                   for node in G.nodes()]

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=NODE_SIZE * area,
        alpha=0.9,
        edgecolors='black',
        linewidths=0.5 * mark,  # Half width border
        ax=ax
    )

    # Draw edges
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

        # Draw edge labels for strong connections
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

    # Draw node labels
    labels = {node: get_short_label(node) for node in G.nodes()}
    nx.draw_networkx_labels(
        G, pos,
        labels=labels,
        font_size=NODE_LABEL_PT,
        font_weight='bold',
        ax=ax
    )

    # Add module title
    ax.set_title(MODULE_DISPLAY[module_id], pad=3 * mark)

    # Set axis limits to prevent label clipping
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')


def main():
    """Generate merged network figure."""
    print("=" * 60)
    print("Merged Cell Interaction Networks - Figure 04B")
    print("=" * 60)

    # Get script directory
    script_dir = Path(__file__).parent

    family = style.apply(title_fontsize=7, fontsize_legend=6,
                         legend_fontsize=6)
    MARK = PANEL_H_MM / EARLIER_CANVAS_H_MM
    AREA = MARK ** 2
    print(f"  type set in {family}; body {style.body_pt():g} pt, "
          f"ticks/legend {style.tick_pt():g} pt; MARK {MARK:.4f}")

    # Create figure with 5x1 subplots
    print("\nCreating figure...")
    fig, axes = style.subplots_mm(PANEL_W_MM, PANEL_H_MM, 5, 1)

    # Process each module
    for module_id in range(1, 6):
        print(f"\nProcessing Module {module_id}...")

        # Load Jaccard matrix
        jaccard_path = script_dir / f'Module_{module_id}_jaccard_matrix.csv'
        jaccard_df = pd.read_csv(jaccard_path, index_col=0)
        print(f"  Jaccard matrix: {jaccard_df.shape}")

        # Create network graph
        G = create_network_graph(jaccard_df, SIMILARITY_THRESHOLD)
        print(f"  Network: {len(G.nodes())} nodes, {len(G.edges())} edges")

        # Plot on corresponding axis
        ax = axes[module_id - 1]
        plot_single_network(ax, G, module_id, MARK, AREA)

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
    print(f"\nSaved: {stem}.[svg|pdf|png] at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print("=" * 60)


if __name__ == "__main__":
    main()
