"""
Reviewer 1 point R1.6 - the stratified result drawn on the tissue.

S9B reports the CEACAM-distance association within tertiles of local epithelial
density as three coefficients. Coefficients are not intuitive, so this panel
shows the same thing on one section: for each spot, a line to each of the five
nearest immune-rich spots, which is exactly the quantity the model uses. The
spot itself is excluded from that set, so no distance falls below one 100 um
spot pitch.

Layout is 2 x 3. Rows are CEACAM-high and CEACAM-low spots, columns are the
three density tertiles, so reading across a row gives the trend and reading down
a column gives the contrast at matched density. That is the same comparison S9B
makes, with no statistics in between.

Section choice is by rule, not by eye: of the ten Visium sections, the one whose
within-high-density coefficient is closest to the cohort median is used.

The figure is for the response letter, not for the paper, so it is written to
this module's outputs rather than into the supplementary figure tree.

Inputs : Round_5/02_Preparation_for_Panels/Spatial/CEACAM_Deconvolution/spot_data.csv
Outputs: spatial_distance_map.csv, R_distance_map.[svg|pdf|png]
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import SPATIAL_SPOT_DATA  # noqa: E402

warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

SCALE, CM, DPI = 4, 1 / 2.54, 300
SEED = 0
# Drawing every spot turns the spokes into a wash of colour; a seeded subsample
# keeps individual lines readable and the ink comparable between panels. The
# median printed on each panel is over all of its spots, not the subsample.
MAX_ORIGINS = 90
K_IMMUNE = 5          # the statistic averages the five nearest, so five are drawn
IMMUNE_THRESHOLD = 0.15
SPOT_PITCH_UM = 100.0  # Visium centre-to-centre, used to put the axes in microns

IMMUNE_TYPES = ["B cells", "CD4+ T cells", "CD8+ T cells", "DC cells",
                "Mast cells", "Monocytes/Macrophages", "NK cells",
                "Neutrophils", "Plasma cells"]
TERTILES = ["Low", "Medium", "High"]
GROUPS = ["CEACAM-high", "CEACAM-low"]
# One line colour for both rows: the row is already named, and a dark line on a
# pale background was hard to see next to the red one.
COLOR_LINE = "#B2182B"
COLOR_SPOT = "#000000"        # the spots this panel measures from
COLOR_OTHER = "#e4e4e4"       # every other spot, for the tissue outline
COLOR_IMMUNE = "#2166AC"      # the spots the lines run to

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def representative_section(spot):
    """The section whose within-high-density coefficient is closest to the median."""
    d = spot.copy()
    d["z"] = d.groupby("sample")["distance_to_immune"].transform(
        lambda s: (s - s.mean()) / s.std(ddof=1))
    high = d[d["tertile"] == "High"]
    betas = {}
    for s, g in high.groupby("sample"):
        if len(g) < 30 or g["CEACAM_ratio"].std(ddof=1) == 0:
            continue
        betas[s] = float(smf.ols("z ~ CEACAM_ratio", g).fit().params["CEACAM_ratio"])
    if not betas:
        raise SystemExit("no section has enough high-density spots")
    b = pd.Series(betas)
    return b.sub(b.median()).abs().idxmin(), b


def upright(xy):
    """
    Rotate the section onto its principal axis. Purely presentational: a rigid
    rotation changes no distance, and it stops a diagonal section wasting most
    of the panel.
    """
    c = xy - xy.mean(axis=0)
    v = np.linalg.svd(c, full_matrices=False)[2][0]
    a = np.arctan2(v[1], v[0])
    r = np.array([[np.cos(-a), -np.sin(-a)], [np.sin(-a), np.cos(-a)]])
    return c @ r.T


def microns_per_pixel(xy):
    """Derived from the array geometry rather than a per-slide scale factor."""
    d, _ = cKDTree(xy).query(xy, k=2)
    return SPOT_PITCH_UM / float(np.median(d[:, 1]))


def nearest_immune(tree, xy, k=K_IMMUNE):
    """
    Distances from each spot to its k nearest immune-rich spots, excluding the
    spot itself. Better than half the spots on a section are immune-rich under
    the 0.15 threshold, and each of those is its own nearest neighbour at zero,
    which pulls the mean below one spot pitch - a distance the array geometry
    cannot produce.
    """
    d, idx = tree.query(xy, k=k + 1)
    self_hit = d[:, 0] < 1e-9
    d = np.where(self_hit[:, None], d[:, 1:], d[:, :-1])
    idx = np.where(self_hit[:, None], idx[:, 1:], idx[:, :-1])
    return d, idx


def main():
    spot = pd.read_csv(SPATIAL_SPOT_DATA).dropna(
        subset=["distance_to_immune", "CEACAM_ratio", "neighborhood_epi_density",
                "sample", "x", "y", "CEACAM_group"])
    spot["tertile"] = pd.qcut(spot["neighborhood_epi_density"], 3, labels=TERTILES)

    section, betas = representative_section(spot)
    s = spot[spot["sample"] == section].copy()
    s["Total_Immune"] = s[IMMUNE_TYPES].sum(axis=1)

    xy = upright(s[["x", "y"]].values)
    s["x"], s["y"] = xy[:, 0], xy[:, 1]
    um = microns_per_pixel(xy)

    immune = s.loc[s["Total_Immune"] > IMMUNE_THRESHOLD, ["x", "y"]].values
    if len(immune) < K_IMMUNE + 1:
        raise SystemExit(f"{section}: only {len(immune)} immune-rich spots")
    tree = cKDTree(immune)
    # The median printed on each panel is the quantity the lines draw, computed
    # here rather than read from the precomputed column, which includes the spot
    # itself at distance zero.
    s["dist_um"] = nearest_immune(tree, xy)[0].mean(axis=1) * um

    rows = []
    for t in TERTILES:
        hi = s[(s["tertile"] == t) & (s["CEACAM_group"] == "CEACAM-high")]["dist_um"]
        lo = s[(s["tertile"] == t) & (s["CEACAM_group"] == "CEACAM-low")]["dist_um"]
        p = stats.mannwhitneyu(hi, lo, alternative="two-sided")[1] if \
            min(len(hi), len(lo)) > 1 else np.nan
        rows.append(dict(section=section, tertile=t,
                         n_ceacam_high=len(hi), n_ceacam_low=len(lo),
                         median_um_ceacam_high=round(float(hi.median()), 1),
                         median_um_ceacam_low=round(float(lo.median()), 1),
                         ratio=round(float(hi.median() / lo.median()), 2),
                         p_two_sided=p))
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "spatial_distance_map.csv", index=False)
    print(f"representative section: {section}  "
          f"(within-high-density beta {betas[section]:.2f}, median {betas.median():.2f})")
    print(table.to_string(index=False))

    _panel(s, immune, tree, table, um, section)


def _panel(s, immune, tree, table, um, section):
    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(2, 3, figsize=(17.5 * SCALE * CM, 7.6 * SCALE * CM))

    x0, x1 = s["x"].min(), s["x"].max()
    y0, y1 = s["y"].min(), s["y"].max()

    for i, grp in enumerate(GROUPS):
        for j, tert in enumerate(TERTILES):
            ax = axes[i, j]
            ax.scatter(s["x"], s["y"], s=2.4 * SCALE, c=COLOR_OTHER,
                       edgecolors="none", zorder=1)
            ax.scatter(immune[:, 0], immune[:, 1], s=2.8 * SCALE, c=COLOR_IMMUNE,
                       alpha=0.45, edgecolors="none", zorder=2)

            sub = s[(s["tertile"] == tert) & (s["CEACAM_group"] == grp)]
            take = sub if len(sub) <= MAX_ORIGINS else sub.iloc[
                np.sort(rng.choice(len(sub), MAX_ORIGINS, replace=False))]
            pts = take[["x", "y"]].values
            _, idx = nearest_immune(tree, pts)
            segs_x, segs_y = [], []
            for p, nb in zip(pts, np.atleast_2d(idx)):
                for k in nb:
                    segs_x += [p[0], immune[k, 0], np.nan]
                    segs_y += [p[1], immune[k, 1], np.nan]
            # Most spokes are one to two spot pitches long once the spot itself
            # is excluded, so they need to be wider than the dots they join.
            ax.plot(segs_x, segs_y, color=COLOR_LINE, linewidth=1.4,
                    alpha=0.8, zorder=4)
            ax.scatter(pts[:, 0], pts[:, 1], s=5.0 * SCALE, c=COLOR_SPOT,
                       edgecolors="white", linewidths=0.3 * SCALE, zorder=5)

            r = table[table["tertile"] == tert].iloc[0]
            med = r["median_um_ceacam_high"] if grp == "CEACAM-high" \
                else r["median_um_ceacam_low"]
            shown = "" if len(take) == len(sub) else f", {len(take)} drawn"
            ax.text(0.02, 0.97, f"median {med:.0f} µm\nn = {len(sub):,}{shown}",
                    transform=ax.transAxes, va="top", ha="left",
                    fontsize=5 * SCALE, color="#333333")

            ax.set_xlim(x0 - 300, x1 + 300)
            ax.set_ylim(y1 + 500, y0 - 300)     # image convention, y downwards
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color("#bbbbbb")
                sp.set_linewidth(0.6)
            if i == 0:
                ax.set_title(f"{tert} epithelial density", fontsize=6.5 * SCALE)
            if j == 0:
                ax.set_ylabel(f"{grp} spots", fontsize=6.5 * SCALE,
                              color="#333333")

    # scale bar, 500 microns, in the lower right panel
    ax = axes[1, 2]
    bar = 500.0 / um
    ax.plot([x1 - bar - 250, x1 - 250], [y1 + 300, y1 + 300],
            color="#333333", linewidth=1.2 * SCALE, solid_capstyle="butt")
    ax.text(x1 - bar / 2 - 250, y1 + 230, "500 µm", ha="center", va="bottom",
            fontsize=5 * SCALE, color="#333333")

    fig.suptitle("Each line joins a spot to one of its five nearest other immune-rich "
                 "spots, the quantity S9B models; one representative section",
                 fontsize=6.5 * SCALE, y=0.985)
    fig.subplots_adjust(left=0.04, right=0.997, top=0.855, bottom=0.015,
                        wspace=0.05, hspace=0.10)
    stem = OUT / "R_distance_map"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
