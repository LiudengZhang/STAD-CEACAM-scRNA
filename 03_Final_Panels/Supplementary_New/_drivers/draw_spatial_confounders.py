#!/usr/bin/env python3
"""Draw Figure S9 from the adopted section-level spatial analysis.

Panels A and B read the section coefficients and their precomputed summaries
from ``spatial_sample_level*.csv``. Panel C redraws the representative map
from the spot table with the exact helpers owned by ``spatial_distance_map``;
its values are checked against ``spatial_distance_map.csv``. No model is fit.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

import _driver_base as base
import panel_style_cns as style

base.apply_style()
SAMPLE = base.analysis("05_R1.6_Spatial_Confounders/scripts/spatial_sample_level.py")
MAP = base.analysis("05_R1.6_Spatial_Confounders/scripts/spatial_distance_map.py")
OUT = base.outputs_of(SAMPLE)
FIG = "S9_Spatial_Confounders"
AB_W, AB_H = 80.0, 52.0
C_W, C_H = 161.0, 68.0
OUTCOMES = list(SAMPLE.OUTCOMES.values())


def frames():
    return dict(
        sections=pd.read_csv(base.require(
            OUT / "spatial_sample_level_sections.csv", "S9 section coefficients")),
        summary=pd.read_csv(base.require(
            OUT / "spatial_sample_level.csv", "S9 section summaries")),
        map_summary=pd.read_csv(base.require(
            OUT / "spatial_distance_map.csv", "S9 map summary")),
    )


def _distance_panel(fr, outcome, panel, stem, save=True):
    # The imported analysis modules set their own matplotlib defaults at
    # module scope. Re-apply the shared publication style after those imports
    # so A/B use the same resolved face and sizes as S8, S10 and S9D/E.
    base.apply_style()
    sec = fr["sections"]
    summ = fr["summary"]
    selected = sec[sec["outcome"] == outcome]
    means = summ[summ["outcome"] == outcome].set_index("stratum")
    fig, ax = style.subplots_mm(AB_W, AB_H)
    rng = np.random.default_rng(0)
    for i, tert in enumerate(SAMPLE.TERTILES):
        values = selected.loc[selected["stratum"] == tert,
                              "ceacam_coef"].to_numpy()
        jitter = rng.uniform(-0.11, 0.11, len(values))
        ax.scatter(np.full(len(values), i) + jitter, values,
                   s=(1.25 * style.PT_PER_MM) ** 2, color=SAMPLE.COLOR,
                   alpha=0.75, edgecolors="white",
                   linewidths=style.EDGE_PT, zorder=3)
        row = means.loc[tert]
        ax.plot([i - 0.25, i + 0.25], [row["mean"]] * 2,
                color="#333333", linewidth=style.RULE_PT, zorder=4)
        ax.plot([i, i], [row["ci_lo"], row["ci_hi"]],
                color="#333333", linewidth=style.RULE_PT, zorder=4)
    ax.axhline(0, color="#999999", linewidth=style.RULE_PT, zorder=1)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.16 * (hi - lo))
    for i, tert in enumerate(SAMPLE.TERTILES):
        row = means.loc[tert]
        label, kw = style.p_text_kw(float(row["p_wilcoxon"]))
        ax.text(i, 0.985, label,
                transform=ax.get_xaxis_transform(), ha="center", va="top", **kw)
    ax.set_xticks(range(3), ["Low", "Medium", "High"])
    ax.set_xlabel("Local epithelial density")
    ax.set_ylabel("CEACAM ratio coefficient, per section")
    ax.set_title(outcome)
    ax.set_xlim(-0.55, 2.55)
    ax.tick_params(axis="both", width=style.RULE_PT, length=2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    base.fit(fig)
    if save:
        base.save(fig, FIG, panel, stem)
    return fig


def draw_A(fr, save=True):
    return _distance_panel(fr, OUTCOMES[0], "S9_A",
                           "S9_A_immune_distance_by_density", save)


def draw_B(fr, save=True):
    return _distance_panel(fr, OUTCOMES[1], "S9_B",
                           "S9_B_stroma_distance_by_density", save)


def _map_frame(expected):
    spot = pd.read_csv(base.require(MAP.SPATIAL_SPOT_DATA, "spatial spot table")).dropna(
        subset=["distance_to_immune", "CEACAM_ratio", "neighborhood_epi_density",
                "sample", "x", "y", "CEACAM_group"])
    spot["tertile"] = pd.qcut(spot["neighborhood_epi_density"], 3,
                              labels=MAP.TERTILES)
    section = expected["section"].iloc[0]
    section_spots = spot[spot["sample"] == section].copy()
    section_spots["Total_Immune"] = section_spots[MAP.IMMUNE_TYPES].sum(axis=1)
    xy = MAP.upright(section_spots[["x", "y"]].to_numpy())
    section_spots["x"], section_spots["y"] = xy[:, 0], xy[:, 1]
    um = MAP.microns_per_pixel(xy)
    immune = section_spots.loc[
        section_spots["Total_Immune"] > MAP.IMMUNE_THRESHOLD, ["x", "y"]].to_numpy()
    tree = cKDTree(immune)
    section_spots["dist_um"] = MAP.nearest_immune(tree, xy)[0].mean(axis=1) * um
    observed = []
    for tert in MAP.TERTILES:
        high = section_spots[(section_spots["tertile"] == tert) &
                             (section_spots["CEACAM_group"] == "CEACAM-high")]["dist_um"]
        low = section_spots[(section_spots["tertile"] == tert) &
                            (section_spots["CEACAM_group"] == "CEACAM-low")]["dist_um"]
        observed.append(dict(
            section=section, tertile=tert, n_ceacam_high=len(high),
            n_ceacam_low=len(low), median_um_ceacam_high=round(high.median(), 1),
            median_um_ceacam_low=round(low.median(), 1),
            ratio=round(high.median() / low.median(), 2),
            p_two_sided=stats.mannwhitneyu(high, low, alternative="two-sided")[1]))
    observed = pd.DataFrame(observed)
    columns = [c for c in expected.columns if c != "section"]
    pd.testing.assert_frame_equal(observed[columns], expected[columns],
                                  check_exact=False, rtol=1e-12, atol=1e-12)
    return section_spots, immune, tree, um


def draw_C(fr, save=True):
    base.apply_style()
    expected = fr["map_summary"]
    section_spots, immune, tree, um = _map_frame(expected)
    rng = np.random.default_rng(MAP.SEED)
    # Each map has a separate header axis.  The earlier version placed the
    # median and spot-count strings inside the maps, where they crossed spots
    # and nearest-neighbour lines.  The four-row grid makes that separation
    # structural: header, map, header, map.  Nothing can be drawn over the
    # annotation because the axes do not share a box.
    fig = style.figure_mm(C_W, C_H)
    gs = fig.add_gridspec(
        4, 3, height_ratios=[0.25, 1.0, 0.17, 1.0],
        hspace=0.05, wspace=0.05)
    headers = np.empty((2, 3), dtype=object)
    axes = np.empty((2, 3), dtype=object)
    for j in range(3):
        headers[0, j] = fig.add_subplot(gs[0, j])
        axes[0, j] = fig.add_subplot(gs[1, j])
        headers[1, j] = fig.add_subplot(gs[2, j])
        axes[1, j] = fig.add_subplot(gs[3, j])
        headers[0, j].set_axis_off()
        headers[1, j].set_axis_off()
    x0, x1 = section_spots["x"].min(), section_spots["x"].max()
    y0, y1 = section_spots["y"].min(), section_spots["y"].max()
    for i, group in enumerate(MAP.GROUPS):
        for j, tert in enumerate(MAP.TERTILES):
            ax = axes[i, j]
            ax.scatter(section_spots["x"], section_spots["y"], s=1.0,
                       color=MAP.COLOR_OTHER, edgecolors="none", zorder=1)
            ax.scatter(immune[:, 0], immune[:, 1], s=1.2,
                       color=MAP.COLOR_IMMUNE, alpha=0.45,
                       edgecolors="none", zorder=2)
            subset = section_spots[(section_spots["tertile"] == tert) &
                                   (section_spots["CEACAM_group"] == group)]
            take = subset if len(subset) <= MAP.MAX_ORIGINS else subset.iloc[
                np.sort(rng.choice(len(subset), MAP.MAX_ORIGINS, replace=False))]
            points = take[["x", "y"]].to_numpy()
            _, indices = MAP.nearest_immune(tree, points)
            sx, sy = [], []
            for point, neighbours in zip(points, np.atleast_2d(indices)):
                for k in neighbours:
                    sx.extend([point[0], immune[k, 0], np.nan])
                    sy.extend([point[1], immune[k, 1], np.nan])
            ax.plot(sx, sy, color=MAP.COLOR_LINE, linewidth=style.RULE_PT,
                    alpha=0.8, zorder=4)
            ax.scatter(points[:, 0], points[:, 1], s=3.0,
                       color=MAP.COLOR_SPOT, edgecolors="white",
                       linewidths=style.EDGE_PT, zorder=5)
            row = expected[expected["tertile"] == tert].iloc[0]
            median = row["median_um_ceacam_high"] if group == "CEACAM-high" \
                else row["median_um_ceacam_low"]
            shown = "" if len(take) == len(subset) else f"; {len(take)} drawn"
            header = headers[i, j]
            header.text(
                0.02, 0.12,
                f"median {median:.0f} µm; n = {len(subset):,}{shown}",
                transform=header.transAxes, va="bottom", ha="left",
                fontsize=style.tick_pt(), color="#333333")
            ax.set_xlim(x0 - 300, x1 + 300)
            ax.set_ylim(y1 + 500, y0 - 300)
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("#bbbbbb")
                spine.set_linewidth(style.RULE_PT)
            if i == 0:
                header.text(
                    0.5, 0.98, f"{tert} epithelial density",
                    transform=header.transAxes, va="top", ha="center",
                    fontsize=style.body_pt(), color="#333333")
            if j == 0:
                ax.set_ylabel(f"{group} spots", color="#333333")
    ax = axes[1, 2]
    bar = 500.0 / um
    ax.plot([x1 - bar - 250, x1 - 250], [y1 + 300, y1 + 300],
            color="#333333", linewidth=style.RULE_PT, solid_capstyle="butt")
    ax.text(x1 - bar / 2 - 250, y1 + 230, "500 µm", ha="center", va="bottom",
            fontsize=style.tick_pt(), color="#333333")
    base.fit(fig)
    if save:
        base.save(fig, FIG, "S9_C", "S9_C_representative_spatial_map")
    return fig


def _check(fr):
    assert set(fr["sections"]["outcome"]) == set(OUTCOMES)
    assert set(fr["summary"]["outcome"]) == set(OUTCOMES)
    _map_frame(fr["map_summary"])
    print("S9: section summaries loaded; representative-map values match "
          "spatial_distance_map.csv")
    return 0


def main():
    fr = frames()
    if "--check" in sys.argv[1:]:
        return _check(fr)
    only = [a for a in sys.argv[1:] if not a.startswith("-")]
    panels = {"S9_A": draw_A, "S9_B": draw_B, "S9_C": draw_C}
    for label, draw in panels.items():
        if not only or label in only:
            draw(fr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
