#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""08_pathway_dotplot.py - Pathway enrichment dotplot"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp
import yaml
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent
CONFIG_PATH = BASE_DIR / "00_Config" / "config.yaml"
SIG_DIR = BASE_DIR / "03_Output" / "03_MetaPrograms" / "metaprogram_signatures"
OUTPUT_DIR = BASE_DIR / "03_Output" / "05_Figures"

with open(CONFIG_PATH, 'r') as f:
    config = yaml.safe_load(f)

gene_sets = config['enrichment']['gene_sets']
top_n = config['enrichment']['top_n_pathways']
padj_cutoff = config['enrichment']['padj_cutoff']

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

sig_files = list(SIG_DIR.glob("MP*_genes.txt"))
all_enrichments = []

for sig_file in sorted(sig_files):
    mp_id = sig_file.stem.replace("_genes", "")
    with open(sig_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]

    for gene_set_lib in gene_sets:
        try:
            enr = gp.enrichr(gene_list=genes, gene_sets=gene_set_lib, organism='Human', outdir=None, cutoff=1.0)
            results = enr.results
            if results is not None and len(results) > 0:
                results['metaprogram'] = mp_id
                results['gene_set_library'] = gene_set_lib
                sig_results = results[results['Adjusted P-value'] < padj_cutoff]
                if len(sig_results) > 0:
                    all_enrichments.append(sig_results.head(top_n))
        except Exception:
            continue

if len(all_enrichments) == 0:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, 'No significant pathway enrichments found\n(Adjusted P-value < 0.05)', ha='center', va='center', fontsize=14)
    ax.axis('off')
    plt.savefig(OUTPUT_DIR / "08_pathway_enrichment_dotplot.png", dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
else:
    combined_df = pd.concat(all_enrichments, ignore_index=True)
    combined_df.to_csv(OUTPUT_DIR / "08_pathway_enrichment_results.csv", index=False)

    top_per_mp = []
    for mp in combined_df['metaprogram'].unique():
        mp_results = combined_df[combined_df['metaprogram'] == mp].sort_values('Adjusted P-value')
        top_per_mp.append(mp_results.head(5))
    plot_df = pd.concat(top_per_mp, ignore_index=True)

    def shorten_name(name, max_len=50):
        return name[:max_len-3] + '...' if len(name) > max_len else name

    plot_df['Term_short'] = plot_df['Term'].apply(shorten_name)
    pivot_pval = plot_df.pivot_table(values='Adjusted P-value', index='Term_short', columns='metaprogram', aggfunc='min')
    pivot_odds = plot_df.pivot_table(values='Odds Ratio', index='Term_short', columns='metaprogram', aggfunc='max')

    fig, ax = plt.subplots(figsize=(max(10, len(pivot_pval.columns) * 1.5), max(8, len(pivot_pval) * 0.4)))

    for i, mp in enumerate(pivot_pval.columns):
        for j, term in enumerate(pivot_pval.index):
            pval = pivot_pval.loc[term, mp]
            odds = pivot_odds.loc[term, mp] if not pd.isna(pivot_odds.loc[term, mp]) else 0
            if pd.isna(pval):
                continue
            size = min(500, max(50, odds * 20))
            color_val = min(10, -np.log10(pval + 1e-10))
            ax.scatter(i, j, s=size, c=color_val, cmap='Reds', vmin=0, vmax=10, edgecolors='black', linewidths=0.5)

    ax.set_xticks(range(len(pivot_pval.columns)))
    ax.set_xticklabels(pivot_pval.columns, rotation=45, ha='right', fontsize=10)
    ax.set_yticks(range(len(pivot_pval.index)))
    ax.set_yticklabels(pivot_pval.index, fontsize=9)
    ax.set_xlabel('Meta-Program', fontsize=12)
    ax.set_ylabel('Pathway', fontsize=12)
    ax.set_title('Pathway Enrichment by Meta-Program', fontsize=14, fontweight='bold')

    sm = plt.cm.ScalarMappable(cmap='Reds', norm=plt.Normalize(vmin=0, vmax=10))
    sm.set_array([])
    plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.02).set_label('-log10(Adj. P-value)', fontsize=10)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "08_pathway_enrichment_dotplot.png", dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

print("Figure 8 complete")
