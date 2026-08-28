#!/usr/bin/env python3
"""
Panel 2E2: Regression diagnostic for C2_Epi_CEACAM6 proportion
  1) Scatter: mean_tumor_score vs C2 proportion with standardized slope (beta)
  2) Boxplot: residuals after regressing out mean_tumor_score, R vs NR
"""
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import EPITHELIAL_H5AD, EPITHELIAL_TUMOR_SCORED_H5AD

# ---------- Data ----------
adata = sc.read_h5ad(EPITHELIAL_H5AD)
adata_ts = sc.read_h5ad(EPITHELIAL_TUMOR_SCORED_H5AD, backed='r')
adata.obs['tumor_score'] = adata_ts.obs['tumor_score'].reindex(adata.obs.index)

pre_mask = adata.obs['Treatment phase'] == 'Pre'
valid_mask = adata.obs['stomach_pre_grouping'].isin(['Responsed', 'No-response'])
sub = adata[pre_mask & valid_mask].copy()
sub.obs['is_C2'] = (sub.obs['minor_cell_state'] == 'C2_Epi_CEACAM6').astype(int)

df = sub.obs.groupby('sample', observed=True).agg(
    C2_proportion=('is_C2', lambda x: x.mean() * 100),
    mean_tumor_score=('tumor_score', 'mean'),
    group=('stomach_pre_grouping', 'first'),
).reset_index().dropna()

x = df['mean_tumor_score'].values
y = df['C2_proportion'].values

slope, intercept, r, p_fit, _ = stats.linregress(x, y)
beta = r                  # standardized slope == pearson r for simple OLS
r_squared = r ** 2
df['residual'] = y - (slope * x + intercept)

R_res = df[df['group'] == 'Responsed']['residual'].values
NR_res = df[df['group'] == 'No-response']['residual'].values
stat_r, p_res = stats.mannwhitneyu(NR_res, R_res, alternative='two-sided')

print(f"beta (standardized slope) = {beta:.4f}")
print(f"R^2 = {r_squared:.4f}   raw slope = {slope:.3f}   p_fit = {p_fit:.4f}")
print(f"Residual boxplot: R mean={R_res.mean():.2f}  NR mean={NR_res.mean():.2f}  p={p_res:.4f}")

# ---------- Plotting ----------
SCALE = 4
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7 * SCALE,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
CM_TO_INCH = 1 / 2.54
COLORS = {'Responsed': '#0072B2', 'No-response': '#D55E00'}
MEDIAN_COLORS = {'Responsed': '#005689', 'No-response': '#A34700'}
OUT = Path(__file__).parent

# ===== Figure 1: scatter with beta =====
fig, ax = plt.subplots(figsize=(4.0 * SCALE * CM_TO_INCH, 3.5 * SCALE * CM_TO_INCH))

for grp, label in [('Responsed', 'R'), ('No-response', 'NR')]:
    sel = df['group'] == grp
    ax.scatter(df.loc[sel, 'mean_tumor_score'], df.loc[sel, 'C2_proportion'],
               s=40 * SCALE, c=COLORS[grp], edgecolors='white',
               linewidths=0.4 * SCALE, label=label, zorder=3)

xs = np.linspace(x.min(), x.max(), 100)
ax.plot(xs, slope * xs + intercept, color='black',
        linewidth=0.8 * SCALE, linestyle='--', zorder=2)

ax.text(0.03, 0.97,
        f'$\\beta$ = {beta:.3f}\n$R^2$ = {r_squared:.3f}\np = {p_fit:.3f}',
        transform=ax.transAxes, ha='left', va='top', fontsize=5 * SCALE)

ax.set_xlabel('Mean tumor score', fontsize=6 * SCALE)
ax.set_ylabel('CEACAM5/6 Epi (%)', fontsize=6 * SCALE)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.5 * SCALE)
ax.spines['bottom'].set_linewidth(0.5 * SCALE)
ax.tick_params(axis='both', labelsize=5 * SCALE,
               width=0.5 * SCALE, length=3 * SCALE)
ax.legend(fontsize=5 * SCALE, frameon=False, loc='upper right')
fig.subplots_adjust(left=0.22, right=0.96, top=0.95, bottom=0.20)
for ext in ('png', 'svg', 'pdf'):
    plt.savefig(OUT / f'tumor_score_vs_c2_scatter.{ext}',
                dpi=300, facecolor='white')
plt.close()

# ===== Figure 2: residual boxplot =====
fig, ax = plt.subplots(figsize=(3.2 * SCALE * CM_TO_INCH, 3.0 * SCALE * CM_TO_INCH))

data = [R_res, NR_res]
bp = ax.boxplot(data, positions=[1, 2], widths=0.6, patch_artist=True,
                boxprops=dict(linewidth=0.5),
                whiskerprops=dict(color='black', linewidth=0.5),
                capprops=dict(color='black', linewidth=0.5),
                flierprops=dict(marker='o', markerfacecolor='white',
                                markersize=4, markeredgecolor='black',
                                markeredgewidth=0.5))
bp['boxes'][0].set_facecolor(COLORS['Responsed'])
bp['boxes'][1].set_facecolor(COLORS['No-response'])
bp['medians'][0].set_color(MEDIAN_COLORS['Responsed'])
bp['medians'][0].set_linewidth(0.8)
bp['medians'][1].set_color(MEDIAN_COLORS['No-response'])
bp['medians'][1].set_linewidth(0.8)

rng = np.random.default_rng(0)
for i, vals in enumerate(data, start=1):
    jitter = rng.uniform(-0.08, 0.08, size=len(vals))
    ax.scatter(np.full_like(vals, i) + jitter, vals,
               s=15 * SCALE, color='black', zorder=4)

y_max = max(np.max(R_res), np.max(NR_res))
y_min = min(np.min(R_res), np.min(NR_res))
span = y_max - y_min
y_bracket = y_max + span * 0.12
ax.plot([1, 1, 2, 2],
        [y_bracket, y_bracket + span * 0.04, y_bracket + span * 0.04, y_bracket],
        'k-', linewidth=0.5)
ax.text(1.5, y_bracket + span * 0.06,
        f'P = {p_res:.3f}',
        ha='center', va='bottom', fontsize=6 * SCALE)

ax.set_title('CEACAM5/6 Epithelial\n(tumor-score adjusted)',
             fontsize=7 * SCALE, fontweight='normal')
ax.set_ylabel('Residual proportion (%)', fontsize=6 * SCALE)
ax.set_xticks([1, 2])
ax.set_xticklabels(['R', 'NR'], fontsize=6 * SCALE)
ax.set_ylim(y_min - span * 0.15, y_max + span * 0.35)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.tick_params(axis='both', labelsize=6 * SCALE, width=0.5, length=4)
fig.subplots_adjust(left=0.28, right=0.95, top=0.85, bottom=0.15)
for ext in ('png', 'svg', 'pdf'):
    plt.savefig(OUT / f'c2_proportion_pre_boxplot_adjusted.{ext}',
                dpi=300, facecolor='white', bbox_inches='tight')
plt.close()

print(f"\nSaved to {OUT}")
