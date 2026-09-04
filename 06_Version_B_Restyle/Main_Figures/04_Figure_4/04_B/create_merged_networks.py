#!/usr/bin/env python3
"""
Figure 4 panel B, RESTYLED (Version B) - the five module networks, one per row,
nodes laid out on a circle and edges weighted by Jaccard similarity.

Version A is
`03_Revised_Panels/Main_Figures/04_Figure_4/04_B/create_merged_networks.py`
and is frozen. This is a copy of it with only the changes PANEL_SPEC.md allows:
the type comes from `00_Config/panel_style_cns.py`, the canvas is the
millimetre box the panel prints in, non-type point sizes are rescaled by MARK,
margins are millimetres, and the save is `style.save_panel`.

Every Jaccard matrix, every threshold, every node selection, every layout
position and every label string is Version A's. The drawing code is the same
code.

  printed panel  Figure 4 B     (PROVENANCE.csv; NOT inferred from "04_B")
  printed rect   21.8 x 94.3 mm    (panel_rects.csv)
  Version B box  50.0 x 178.0 mm

DETERMINISM
    Checked before anything was changed, because this panel has had a
    layout-moves-between-runs problem in the past. It does not now, and it
    cannot: `nx.circular_layout` places nodes by their order in the graph, the
    order comes from the CSV index, and the `MAX_NODES` trim uses `sorted(...)`,
    which is stable. There is no RNG in the drawing path and therefore no seed
    to set. `compare_panel_content.py` run with the Version A script on both
    sides prints CONTENT IDENTICAL over all five axes.

MARK - and why SMALL_PT is 4 and not 2
    Version A drew 7 x 12 cm with no SCALE constant and set three type sizes:
    `NODE_FONT_SIZE = 2` on the node labels, `TITLE_FONT_SIZE = 4` on the module
    titles, and `EDGE_FONT_SIZE = 5` on edge labels that are never drawn
    (`EDGE_LABEL_THRESHOLD = 1.1` and a Jaccard cannot exceed 1.0). The
    assembler cropped the figure to its ink and fitted the resulting
    63.6 x 333.0 pt SVG into 21.8 x 94.3 mm, a fit of 0.8027: the module titles
    printed at 3.21 pt and the node labels at 1.61 pt.

    Read literally, PANEL_SPEC's "smallest body type" is the 2 pt node label and
    MARK would be 7 / 2 = 3.5. That is degenerate here. It would set the node
    discs, whose area is `NODE_SIZE = 150`, at 3.5^2 x 150 = 1837 pt^2, i.e.
    15.1 mm across, on a network circle about 30 mm wide - eight of them on a
    circle of 76 mm circumference, needing 121 mm. The panel would be unreadable
    in a way Version A is not.

    The reason is that the derivation assumes the panel box grows by the same
    ratio the type does. Every other Figure 4 panel can do that. This one cannot:
    at 4.36x the printed box it would be 95 x 411 mm, taller than any page. So
    the reference is the panel's smallest *legible* body type, the 3.21 pt module
    title, whose nominal value is 4:

        SMALL_PT = 4, SCALE = 1
        MARK = tick_pt / (SMALL_PT * SCALE) = 7 / 4 = 1.75
        AREA = MARK ** 2 = 3.0625

    which puts the node disc at 7.6 mm on a ~30 mm network, against 3.47 mm on
    the ~18.9 mm network as printed - the same order of proportion, 0.25 against
    0.18. This is the one judgement call in Figure 4 that departs from the
    letter of PANEL_SPEC, and it is recorded here rather than buried.

    WHAT IS LOST, STATED PLAINLY. In the published panel each node label sits
    *inside* its disc: 2 pt type on a 12.25 pt disc. Here the label is 7 pt and
    the disc 21.4 pt, so the label is wider than the disc it names and overhangs
    it. Restoring the published relation needs a disc of at least 33 pt, which
    is MARK 2.7 or more, and eight of those need a network 67 mm across - a
    95 x 411 mm panel. There is no size on any page at which this panel both
    reaches 7 pt and keeps its printed proportions. Legibility was chosen over
    proportion. Nothing plotted moved; the content gate passes.

THE BOX - a very tall, very narrow panel that has to get wider
    21.8 mm cannot carry a node label at 7 pt: 'C3_Prolif' is about 11.6 mm
    wide, over half the published width. The label geometry sets both
    dimensions. Two labels on adjacent nodes of an eight-node circle are
    separated vertically by 0.293 r, so the network square must be at least
    ~24 mm across for them to clear each other; the labels then run past the
    axes' own box, which is why the figure is 50 mm wide while the network
    squares are ~30 mm - the extra 10 mm a side is the room the labels need,
    and it keeps every one of them inside the canvas (`overflow_mm` returns
    zeros). Five stacked squares plus their titles is 178 mm. The published
    aspect (0.23) becomes 0.28: preserving it would put the labels back on top
    of each other, which the brief explicitly releases this panel from.

Input : this directory / Module_{1..5}_jaccard_matrix.csv
        Version A reads these from beside its own script and the path is part of
        the code PANEL_SPEC forbids changing, so the five files (and the unused
        module_mappings.csv that sits with them) were copied here byte for byte.
        Nothing was moved or deleted.
Output: this directory / merged_networks_k5.{svg,pdf,png}
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "00_Config"))
import panel_style_cns as style              # noqa: E402

# Version A's canvas convention, kept only so MARK can be derived from it.
SCALE = 1
SMALL_PT = 4.0                # Version A's `TITLE_FONT_SIZE`; see the docstring

PRINTED_MM = (21.8, 94.3)     # published rect, panel_rects.csv
PANEL_W_MM, PANEL_H_MM = 50.0, 178.0
MARGIN = dict(left=1.0, right=1.0, top=4.5, bottom=1.0)
HSPACE = 0.30                 # between the five rows, as a fraction of a row

# Parameters
DPI = 600
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


def plot_single_network(ax, G, module_id, MARK, AREA):
    """Plot a single network on a given axis."""
    # Limit nodes if needed
    if len(G.nodes()) > MAX_NODES:
        degrees = dict(G.degree())
        top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:MAX_NODES]
        G = G.subgraph(top_nodes).copy()

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
        node_size=NODE_SIZE * AREA,
        alpha=0.9,
        edgecolors='black',
        linewidths=0.5 * MARK,  # Half width border
        ax=ax
    )

    # Draw edges
    edges = G.edges()
    if len(edges) > 0:
        weights = [G[u][v]['weight'] for u, v in edges]
        edge_widths = [w * EDGE_WIDTH_MULTIPLIER * MARK for w in weights]

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
        font_size=style.tick_pt(),
        ax=ax
    )

    # Add module title
    ax.set_title(MODULE_DISPLAY[module_id], pad=3 * MARK)

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

    family = style.apply()
    MARK = style.tick_pt() / (SMALL_PT * SCALE)
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

    style.margins_mm(fig, **MARGIN, hspace=HSPACE)
    over = style.overflow_mm(fig)
    if max(over) > 0.05:
        print(f"  WARNING ink outside the canvas (l,r,b,t mm): "
              f"{tuple(round(v, 2) for v in over)}")
    style.save_panel(fig, script_dir / 'merged_networks_k5')
    print(f"\nSaved: {script_dir / 'merged_networks_k5'}.[svg|pdf|png] "
          f"at {PANEL_W_MM} x {PANEL_H_MM} mm")
    print("=" * 60)


if __name__ == "__main__":
    main()
