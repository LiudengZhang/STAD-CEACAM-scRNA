"""Shared drawing code for Figure 3's spatial panels.

    spot_map        F, G, K, L, M: one Visium section, spots coloured by one
                    column of spot_data.csv, with a colour bar; all five at
                    one scale
    paired_boxplot  H, I, J: per-section means in CEACAM5/6-low and -high
                    spots, paired within section

Each printed panel keeps its own script, which reads the spot data and calls
one of these with its column, colours, title and output name.
"""

import matplotlib.pyplot as plt
import pandas as pd

from . import style
from .layout import pin_frame_mm
from .boxes import box_xlim, printed_and_two_sided

# Spot marker area as a multiple of the panel's AREA; small enough that
# neighbouring spots stay separate.
SPOT_AREA = 2.5

# Paired boxplot colours: CEACAM5/6-low, then -high.
LOW_COLOR = '#1f77b4'
HIGH_COLOR = '#d62728'


# Map geometry, mm, shared by the five maps so they print at one scale: an
# equal-aspect frame MAP_TOP_MM below the slot top and MAP_LEFT_MM right of the
# slot named by `under` (K under F, L under G, M under H), with the title flush
# left above it.
MAP_FRAME_MM = (18.0, 24.0)
MAP_TOP_MM = 5.0
MAP_LEFT_MM = 3.4             # the letter cell is 2.87 mm wide
MAP_CBAR_GAP_MM = 0.6
MAP_CBAR_W_MM = 1.4


def spot_map(panel, sample_df, column, cmap, title, stem, under=None,
             cbar_label=None):
    """One section's spots, coloured by `column`, with a colour bar.

    The frame is MAP_FRAME_MM at MAP_TOP_MM, MAP_LEFT_MM right of the left
    edge of the slot of printed panel `under` (default: its own), so maps in
    different rows share left edges. Colour limits are the data's own. Neither
    axis carries a label or a ruler: the coordinates are Visium pixel
    positions within one section.
    """
    import slots
    fw, fh = MAP_FRAME_MM
    fig = panel.canvas()
    W, H = panel.w_mm, panel.h_mm
    own_x0 = slots.rect_mm(panel.figure, panel.letter)[0]
    ref_x0 = slots.rect_mm(panel.figure, under or panel.letter)[0]
    left = ref_x0 + MAP_LEFT_MM - own_x0
    if left < panel.letter_cell[0] + 0.3:
        raise RuntimeError(f"{panel.letter}: the frame at {left:.2f} mm would "
                           f"put the title into the letter cell")
    top = MAP_TOP_MM
    ax = fig.add_axes([left / W, 1 - (top + fh) / H, fw / W, fh / H])

    sc = ax.scatter(sample_df['x'], sample_df['y'],
                    c=sample_df[column],
                    cmap=cmap, s=SPOT_AREA * panel.area, alpha=0.8)
    # Flush-left title, set through the default title so ax.get_title() finds it.
    t = ax.set_title(title, fontweight='normal')
    t.set_horizontalalignment('left')
    t.set_x(0.0)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(style.RULE_PT)
        spine.set_visible(True)
    # Equal aspect within a fixed frame: one scale for every map.
    ax.set_aspect('equal', adjustable='datalim')

    cax = fig.add_axes([(left + fw + MAP_CBAR_GAP_MM) / W,
                        1 - (top + fh) / H, MAP_CBAR_W_MM / W, fh / H])
    cbar = fig.colorbar(sc, cax=cax)
    cbar.solids.set_rasterized(False)
    cbar.ax.tick_params(width=style.RULE_PT, length=4 * panel.mark)
    cbar.outline.set_linewidth(style.RULE_PT)
    if cbar_label:
        cbar.set_label(cbar_label, fontsize=style.tick_pt(), labelpad=1.0)
    return panel.finish(fig, stem)


def section_means(df, col):
    """One mean of `col` per section per CEACAM5/6 group, paired sections only.

    A section that contributed spots to only one of the two groups is dropped,
    because the test printed on the panel is paired within section.
    """
    rows = []
    for sample in df['sample'].unique():
        sample_df = df[df['sample'] == sample]
        for group in ['CEACAM-high', 'CEACAM-low']:
            group_df = sample_df[sample_df['CEACAM_group'] == group]
            if len(group_df) > 0:
                rows.append({'sample': sample, 'group': group,
                             col: group_df[col].mean()})
    out = pd.DataFrame(rows)
    paired = out.groupby('sample').filter(lambda x: len(x) == 2)['sample'].unique()
    out = out[out['sample'].isin(paired)]
    return out.sort_values(['sample', 'group'])


def paired_boxplot(panel, sample_data, col, title, ylabel, stem,
                   y_max_limit=None, tight_y=False, fit_kw=None, exact=False):
    """CEACAM5/6-low against -high, one line per section joining its own pair.

    The panel prints the one-tailed Wilcoxon signed-rank P (high > low); the
    legend reports the two-sided value.

    `tight_y` shortens the y ticks and their pad. `fit_kw` is passed to
    `panel.fit`.
    """
    fig, ax = panel.subplots()

    low_data = sample_data[sample_data['group'] == 'CEACAM-low'][col].values
    high_data = sample_data[sample_data['group'] == 'CEACAM-high'][col].values

    pval, pval_two_sided = printed_and_two_sided(high_data, low_data,
                                                 kind="wilcoxon")

    bp = ax.boxplot([low_data, high_data], positions=[1, 2], widths=0.5,
                    patch_artist=True,
                    boxprops=dict(linewidth=style.RULE_PT),
                    whiskerprops=dict(linewidth=style.RULE_PT),
                    capprops=dict(linewidth=style.RULE_PT),
                    medianprops=dict(linewidth=style.RULE_PT),
                    flierprops=dict(markersize=6.0 * panel.mark,
                                    markeredgewidth=style.EDGE_PT))
    bp['boxes'][0].set_facecolor(LOW_COLOR)
    bp['boxes'][1].set_facecolor(HIGH_COLOR)
    for box in bp['boxes']:
        box.set_alpha(0.7)
        box.set_edgecolor('black')
    for median in bp['medians']:
        median.set_color('black')

    # The pairing is the point of the panel, so the per-section points stay.
    for i in range(len(low_data)):
        ax.plot([1, 2], [low_data[i], high_data[i]], 'k-', alpha=0.3,
                linewidth=style.RULE_PT)

    ax.scatter([1] * len(low_data), low_data, color=LOW_COLOR,
               s=40 * panel.area, zorder=3,
               edgecolor='black', linewidth=style.EDGE_PT)
    ax.scatter([2] * len(high_data), high_data, color=HIGH_COLOR,
               s=40 * panel.area, zorder=3,
               edgecolor='black', linewidth=style.EDGE_PT)

    data_y_max = max(max(low_data), max(high_data))
    y_range = data_y_max - min(min(low_data), min(high_data))
    bracket_y = data_y_max + 0.08 * y_range

    ax.plot([1, 1, 2, 2],
            [bracket_y, bracket_y + 0.02 * y_range,
             bracket_y + 0.02 * y_range, bracket_y],
            color='black', linewidth=style.RULE_PT)

    # A star at STAR_PT, a value at the tick size.
    pval_text, p_kw = style.p_text_kw(pval, exact=exact)
    ax.text(1.5, bracket_y + 0.04 * y_range, pval_text,
            ha='center', va='bottom', **p_kw)

    # One extra point of pad keeps the title off the top tick label.
    ax.set_title(title, pad=plt.rcParams['axes.titlepad'] + 1.0)
    ax.set_ylabel(ylabel)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Low', 'High'])
    ax.set_xlabel('CEACAM5/6', labelpad=1)
    # 0.3 box widths of paper each side, so neither box stands on the y axis.
    ax.set_xlim(*box_xlim([1, 2], 0.5, clear=0.3))
    ax.tick_params(axis='both', width=style.RULE_PT, length=4 * panel.mark)
    # Keeps the lowest y tick label clear of the left group name.
    ax.tick_params(axis='x', pad=2.5)
    if tight_y:
        ax.tick_params(axis='y', pad=0.6, length=1.5)

    for spine in ax.spines.values():
        spine.set_linewidth(style.RULE_PT)

    if y_max_limit is not None:
        ax.set_ylim(0, y_max_limit)

    panel.fit(fig, **(fit_kw or {}))
    pin_frame_mm(fig, ax, top_mm=5.0, bottom_mm=panel.h_mm - 29.0)
    return panel.finish(fig, stem)
