#!/usr/bin/env python3
"""
Create Merged Cell Interaction Networks (Figure 04 Panel B)

Generates a 1x5 merged network visualization showing all 5 module networks.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path

plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']})

# Parameters - Target: 7cm × 12cm at 300 DPI
FIGURE_WIDTH = 7 / 2.54  # Convert cm to inches
FIGURE_HEIGHT = 12 / 2.54  # Increased height for better spacing
DPI = 600
SIMILARITY_THRESHOLD = 0.01
NODE_SIZE = 150  # Smaller nodes
NODE_FONT_SIZE = 2  # Smaller font
EDGE_WIDTH_MULTIPLIER = 2
EDGE_ALPHA = 0.7
EDGE_LABEL_THRESHOLD = 1.1  # Effectively disables edge labels (max Jaccard is 1.0)
EDGE_FONT_SIZE = 5
TITLE_FONT_SIZE = 4  # For module titles (half size, no bold)

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


def plot_single_network(ax, G, module_id):
    """Plot a single network on a given axis."""
    # Limit nodes if needed
    if len(G.nodes()) > MAX_NODES:
        degrees = dict(G.degree())
        top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:MAX_NODES]
        G = G.subgraph(top_nodes).copy()

    if len(G.nodes()) == 0:
        ax.text(0.5, 0.5, f'{MODULE_DISPLAY[module_id]}\n(No nodes)',
                ha='center', va='center', fontsize=TITLE_FONT_SIZE)
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
        node_size=NODE_SIZE,
        alpha=0.9,
        edgecolors='black',
        linewidths=0.5,  # Half width border
        ax=ax
    )

    # Draw edges
    edges = G.edges()
    if len(edges) > 0:
        weights = [G[u][v]['weight'] for u, v in edges]
        edge_widths = [w * EDGE_WIDTH_MULTIPLIER for w in weights]

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
                font_size=EDGE_FONT_SIZE,
                ax=ax
            )

    # Draw node labels
    labels = {node: get_short_label(node) for node in G.nodes()}
    nx.draw_networkx_labels(
        G, pos,
        labels=labels,
        font_size=NODE_FONT_SIZE,
        font_weight='bold',
        ax=ax
    )

    # Add module title
    ax.set_title(MODULE_DISPLAY[module_id], fontsize=TITLE_FONT_SIZE, pad=3)

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

    # Create figure with 5x1 subplots
    print("\nCreating figure...")
    fig, axes = plt.subplots(5, 1, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))

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
        plot_single_network(ax, G, module_id)

    plt.tight_layout()

    # Save
    output_path = script_dir / 'merged_networks_k5.png'
    fig.savefig(output_path, dpi=DPI, bbox_inches='tight')
    fig.savefig(output_path.with_suffix('.svg'), dpi=DPI, bbox_inches='tight')

    file_size = output_path.stat().st_size / 1024
    print(f"\n✓ Saved: {output_path.name}")
    print(f"  Size: {file_size:.1f} KB")
    print("=" * 60)

    plt.close()


if __name__ == "__main__":
    main()
