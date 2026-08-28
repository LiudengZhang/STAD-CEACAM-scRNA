"""
WP4 / Reviewer 1 point R1.4.

  "The Results section states that MP4 and MP5 are strongly associated with
   CEACAM5/6 but then reports: 'MP4 was significantly depleted in pre-treatment
   non-responders...'. Please clarify this statement, as Figure 2H-I appears to
   show the opposite trend. In addition, the same comparison should be performed
   between pre- and post-treatment responders, as well as between pre- and
   post-treatment non-responders."

Two things are settled here.

1. Direction. The published test compares pre-treatment RESPONDERS against all
   other groups pooled, and MP4 is *lower* in pre-treatment responders. The
   sentence in the manuscript names non-responders, which inverts the result.
   The reviewer is right; the figure is right; the sentence is wrong.

2. The four requested group contrasts, all two-sided.

Input : Round_5/02_Preparation_for_Panels/Metaprogram_Permutation/mp4_permutation_results.json
Output: mp_group_comparisons.csv, mp_direction_report.txt, and panel S8A.
"""

from itertools import combinations
from pathlib import Path
import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import PREPARATION, REVISED_PANELS  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
PANEL_DIR = REVISED_PANELS / "Supplementary_New" / "S8_CEACAM_Metaprogram" / "S8_A"
PANEL_DIR.mkdir(parents=True, exist_ok=True)

SCALE, CM, DPI = 4, 1 / 2.54, 300
GROUPS = ["Pre-R", "Pre-NR", "Post-R", "Post-NR"]
COLORS = {"Pre-R": "#bde0fe", "Pre-NR": "#a2d2ff",
          "Post-R": "#ffcfd2", "Post-NR": "#f1c0e8"}

CONTRASTS = [
    ("Pre-NR", "Pre-R", "Pre-treatment NR vs R (the CEACAM comparison)"),
    ("Post-NR", "Post-R", "Post-treatment NR vs R"),
    ("Post-R", "Pre-R", "Post- vs pre-treatment, responders (requested)"),
    ("Post-NR", "Pre-NR", "Post- vs pre-treatment, non-responders (requested)"),
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def exact_perm_two_sided(a, b):
    """Two-sided exact permutation test on the difference in means."""
    pooled = np.concatenate([a, b])
    obs = a.mean() - b.mean()
    diffs = np.array([
        pooled[list(c)].mean() - pooled[[i for i in range(len(pooled)) if i not in c]].mean()
        for c in combinations(range(len(pooled)), len(a))
    ])
    return obs, float(np.mean(np.abs(diffs) >= abs(obs))), len(diffs)


def main():
    mp = json.loads(
        (PREPARATION / "Metaprogram_Permutation" / "mp4_permutation_results.json").read_text())

    vals = {prog: {g: np.asarray(mp[prog][f"{g}_values"], float) for g in GROUPS}
            for prog in ("S-MP4", "S-MP5")}

    rows = []
    for prog in ("S-MP4", "S-MP5"):
        v = vals[prog]

        # (1) the published contrast, restated in the correct direction
        pre_r = v["Pre-R"]
        others = np.concatenate([v[g] for g in ("Post-R", "Pre-NR", "Post-NR")])
        obs, p2, nperm = exact_perm_two_sided(pre_r, others)
        rows.append(dict(
            program=prog, contrast="Pre-R vs all other groups (published test)",
            group_a="Pre-R", group_b="All others",
            n_a=len(pre_r), n_b=len(others),
            mean_a=pre_r.mean(), mean_b=others.mean(),
            diff=obs, test=f"Exact permutation ({nperm} permutations)",
            p_two_tailed=p2,
            direction=("lower in Pre-R" if obs < 0 else "higher in Pre-R"),
        ))

        # (2) the four pairwise group contrasts, Mann-Whitney two-sided
        for ga, gb, label in CONTRASTS:
            a, b = v[ga], v[gb]
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            r = 2.0 * u / (len(a) * len(b)) - 1.0
            rows.append(dict(
                program=prog, contrast=label, group_a=ga, group_b=gb,
                n_a=len(a), n_b=len(b), mean_a=a.mean(), mean_b=b.mean(),
                diff=a.mean() - b.mean(), test="Mann-Whitney U (two-sided)",
                p_two_tailed=float(p), rank_biserial_r=float(r),
                direction=(f"higher in {ga}" if a.mean() > b.mean() else f"higher in {gb}"),
            ))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "mp_group_comparisons.csv", index=False)

    # ------------------------------------------------------------- report
    L = ["METAPROGRAM DIRECTION AND PRE/POST CONTRASTS - Reviewer 1 point R1.4",
         "=" * 90, ""]
    L.append("PART 1 - THE DIRECTION OF THE PUBLISHED MP4 RESULT")
    L.append("-" * 90)
    for prog in ("S-MP4", "S-MP5"):
        v = vals[prog]
        L.append(f"  {prog} group means:")
        for g in GROUPS:
            L.append(f"     {g:<8} n={len(v[g])}  mean={v[g].mean():.4f}  "
                     f"median={np.median(v[g]):.4f}")
        r = df[(df["program"] == prog) & (df["contrast"].str.startswith("Pre-R vs all"))].iloc[0]
        L.append(f"     Pre-R vs all others: difference = {r['diff']:+.4f} "
                 f"({r['direction']}), two-sided {r['test']} P = {r['p_two_tailed']:.4f}")
        L.append("")
    L.append("  The published test contrasts PRE-TREATMENT RESPONDERS against every")
    L.append("  other group. MP4 is LOWER in pre-treatment responders, i.e. HIGHER in")
    L.append("  non-responders. The manuscript sentence names non-responders and is")
    L.append("  therefore inverted. The reviewer's reading of the figure is correct.")
    L.append("")
    L.append("  Corrected wording:")
    L.append('     "MP4 was significantly depleted in pre-treatment responders relative')
    L.append('      to all other groups (two-sided exact permutation P = '
             f'{df[(df["program"] == "S-MP4") & (df["contrast"].str.startswith("Pre-R vs all"))].iloc[0]["p_two_tailed"]:.3f}), that is,')
    L.append('      enriched in non-responders, consistent with the CEACAM5/6 result."')
    L.append("")
    L.append("")
    L.append("PART 2 - THE REQUESTED GROUP CONTRASTS (all two-sided)")
    L.append("-" * 90)
    for prog in ("S-MP4", "S-MP5"):
        L.append(f"  [{prog}]")
        sub = df[(df["program"] == prog) & (df["test"].str.startswith("Mann-Whitney"))]
        for _, r in sub.iterrows():
            L.append(f"     {r['contrast']}")
            L.append(f"        {r['group_a']} mean={r['mean_a']:.4f} (n={r['n_a']})   "
                     f"{r['group_b']} mean={r['mean_b']:.4f} (n={r['n_b']})")
            L.append(f"        P = {r['p_two_tailed']:.4f}   "
                     f"rank-biserial r = {r['rank_biserial_r']:+.3f}   {r['direction']}")
        L.append("")
    report = "\n".join(L)
    (OUT / "mp_direction_report.txt").write_text(report, encoding="utf-8")
    print(report)

    # -------------------------------------------------------------- panel
    fig, axes = plt.subplots(1, 2, figsize=(9.0 * SCALE * CM, 4.6 * SCALE * CM))
    for ax, prog in zip(axes, ("S-MP4", "S-MP5")):
        v = vals[prog]
        data = [v[g] for g in GROUPS]
        bp = ax.boxplot(data, positions=range(4), widths=0.6, patch_artist=True,
                        showfliers=False,
                        boxprops=dict(linewidth=0.8),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8),
                        medianprops=dict(color="black", linewidth=1.2))
        for patch, g in zip(bp["boxes"], GROUPS):
            patch.set_facecolor(COLORS[g])
            patch.set_edgecolor("#444444")
        rng = np.random.default_rng(0)
        for i, g in enumerate(GROUPS):
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(v[g])), v[g],
                       s=8 * SCALE, c="#333333", zorder=3, alpha=0.8,
                       edgecolors="white", linewidths=0.3 * SCALE)

        # Brackets for the four requested contrasts
        top = max(x.max() for x in data)
        step = 0.075 * top
        y = top + step
        for (ga, gb, _), _ in zip(CONTRASTS, range(4)):
            ia, ib = GROUPS.index(ga), GROUPS.index(gb)
            p = df[(df["program"] == prog) & (df["group_a"] == ga)
                   & (df["group_b"] == gb)]["p_two_tailed"].iloc[0]
            lo, hi = sorted((ia, ib))
            ax.plot([lo, lo, hi, hi], [y, y + step * 0.25, y + step * 0.25, y],
                    color="#444444", linewidth=0.8)
            ax.text((lo + hi) / 2, y + step * 0.3, f"P = {p:.3f}",
                    ha="center", va="bottom", fontsize=4.5 * SCALE)
            y += step * 1.15

        ax.set_xticks(range(4))
        ax.set_xticklabels(GROUPS, fontsize=6 * SCALE)
        ax.set_ylabel(f"{prog} score", fontsize=6 * SCALE)
        ax.set_title(prog, fontsize=7 * SCALE)
        ax.set_ylim(0, y + step)
        ax.tick_params(axis="both", labelsize=5.5 * SCALE, width=0.8, length=3)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.spines["left"].set_linewidth(0.8)
        ax.spines["bottom"].set_linewidth(0.8)

    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.10, wspace=0.28)
    stem = PANEL_DIR / "S8_A_metaprogram_four_groups"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"\nSaved {stem}.[svg|pdf|png]")


if __name__ == "__main__":
    main()
