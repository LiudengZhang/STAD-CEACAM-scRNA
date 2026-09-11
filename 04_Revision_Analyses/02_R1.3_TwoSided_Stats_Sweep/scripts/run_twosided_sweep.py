"""
WP3 / Reviewer 1 point R1.3c.

"Two-sided tests should be used throughout. The authors should demonstrate that
the conclusions remain robust under two-sided testing and report exact P values
rather than only significance thresholds."

This script recomputes every comparison that was originally reported with a
one-tailed test, from the same source data the panel scripts use, and reports
both tails side by side together with an effect size and a bootstrap confidence
interval. The point of the CI is that it does not depend on the choice of tail,
so it is the statistic the revised manuscript leans on.

Effect sizes
  Mann-Whitney U        -> rank-biserial correlation r = 2U/(n1*n2) - 1
  Wilcoxon signed-rank  -> matched-pairs rank-biserial correlation
  exact permutation     -> observed difference in group means

The row for Supplementary Figure S2D (the tumor-content-adjusted CEACAM5/6+
proportion) is computed last, after every other row, so that adding it leaves
each earlier row's bootstrap draws exactly where they were. Order matters here:
the draws come from one shared stream.

Units
  Every row is in the units of its own source, and the label says so where a
  unit exists. The two spatial distance rows are in micrometres as of
  2026-09-10, converted here from the array units of spot_data.csv by the one
  factor 00_Config/spatial_scale.py derives; see the comment at section 5.

Outputs (04_Revision_Analyses/02_R1.3_TwoSided_Stats_Sweep/outputs/)
  twosided_sweep.csv        one row per comparison
  twosided_sweep_report.txt human-readable summary with the verdict per row
"""

from itertools import combinations
from pathlib import Path
import json
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import (  # noqa: E402
    EPITHELIAL_H5AD, EPITHELIAL_DEPOSIT_H5AD, MOMAC_H5AD, FIBROBLAST_H5AD,
    DC_CELLS_H5AD, TCD4_H5AD, TIGER_BAYESPRISM_EPI, TIGER_META, SPATIAL_SPOT_DATA,
    PREPARATION,
)
from spatial_scale import cohort_um_per_unit  # noqa: E402
from shared.expression import expression_source  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)
N_BOOT = 10000
MIN_CELLS = 20

rows = []


# ----------------------------------------------------------------- statistics
def p_floor_mw(n1, n2):
    """
    Smallest two-sided P a Mann-Whitney U test can return at these group sizes.
    With complete separation U is extremal, and the exact null puts mass
    1/C(n1+n2, n1) on each tail. At n = 4 vs 4 that floor is 2/70 = 0.029, so
    the test has almost no resolution below 0.05 regardless of effect size.
    """
    from math import comb
    return 2.0 / comb(n1 + n2, n1)


def p_floor_signrank(n):
    """Smallest two-sided P a Wilcoxon signed-rank test can return at n pairs."""
    return 2.0 / (2.0 ** n)


def hedges_g(a, b, paired=False):
    """
    Standardised mean difference with the small-sample correction. This is the
    common scale for the forest plot, since raw units differ across panels
    (log-normalised expression, staining percent, micrometres, NMF score).
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    if paired:
        d = a - b
        sd = d.std(ddof=1)
        g = d.mean() / sd if sd > 0 else 0.0
        n = len(d)
    else:
        n1, n2 = len(a), len(b)
        sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1))
                     / (n1 + n2 - 2))
        g = (a.mean() - b.mean()) / sp if sp > 0 else 0.0
        n = n1 + n2
    return float(g * (1 - 3 / (4 * n - 9))) if n > 3 else float(g)


def boot_g_ci(a, b, paired=False, n=N_BOOT):
    """Percentile bootstrap CI for Hedges' g."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    out = []
    if paired:
        d = a - b
        idx = RNG.integers(0, len(d), size=(n, len(d)))
        for row in d[idx]:
            out.append(hedges_g(row, np.zeros_like(row), paired=True))
    else:
        ia = RNG.integers(0, len(a), size=(n, len(a)))
        ib = RNG.integers(0, len(b), size=(n, len(b)))
        for xa, xb in zip(a[ia], b[ib]):
            out.append(hedges_g(xa, xb))
    out = np.asarray(out)
    out = out[np.isfinite(out)]
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def boot_ci(a, b, paired=False, n=N_BOOT):
    """Percentile bootstrap CI for mean(a) - mean(b)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if paired:
        d = a - b
        idx = RNG.integers(0, len(d), size=(n, len(d)))
        diffs = d[idx].mean(axis=1)
    else:
        ia = RNG.integers(0, len(a), size=(n, len(a)))
        ib = RNG.integers(0, len(b), size=(n, len(b)))
        diffs = a[ia].mean(axis=1) - b[ib].mean(axis=1)
    return float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def add_mw(label, family, panel, hi, lo, hi_name, lo_name, direction="greater"):
    """Mann-Whitney U, testing hi vs lo. `direction` is the original one-tail."""
    hi, lo = np.asarray(hi, float), np.asarray(lo, float)
    u1, p1 = stats.mannwhitneyu(hi, lo, alternative=direction)
    _, p2 = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    r = 2.0 * u1 / (len(hi) * len(lo)) - 1.0
    ci_lo, ci_hi = boot_ci(hi, lo)
    g_lo, g_hi = boot_g_ci(hi, lo)
    rows.append(dict(
        analysis=label, family=family, panel=panel, test="Mann-Whitney U",
        group_hi=hi_name, group_lo=lo_name, n_hi=len(hi), n_lo=len(lo),
        mean_hi=float(np.mean(hi)), mean_lo=float(np.mean(lo)),
        median_hi=float(np.median(hi)), median_lo=float(np.median(lo)),
        statistic=float(u1), p_one_tailed=float(p1), p_two_tailed=float(p2),
        effect_size_name="rank-biserial r", effect_size=float(r),
        diff_of_means=float(np.mean(hi) - np.mean(lo)),
        ci95_lo=ci_lo, ci95_hi=ci_hi,
        p_two_tailed_floor=p_floor_mw(len(hi), len(lo)),
        hedges_g=hedges_g(hi, lo),
        g_ci95_lo=g_lo, g_ci95_hi=g_hi,
    ))


def add_wilcoxon(label, family, panel, hi, lo, hi_name, lo_name):
    """Wilcoxon signed-rank for the paired spatial comparisons."""
    hi, lo = np.asarray(hi, float), np.asarray(lo, float)
    s1, p1 = stats.wilcoxon(hi, lo, alternative="greater")
    _, p2 = stats.wilcoxon(hi, lo, alternative="two-sided")
    d = hi - lo
    nz = d[d != 0]
    ranks = stats.rankdata(np.abs(nz))
    r = (ranks[nz > 0].sum() - ranks[nz < 0].sum()) / ranks.sum()
    ci_lo, ci_hi = boot_ci(hi, lo, paired=True)
    g_lo, g_hi = boot_g_ci(hi, lo, paired=True)
    rows.append(dict(
        analysis=label, family=family, panel=panel,
        test="Wilcoxon signed-rank (paired)",
        group_hi=hi_name, group_lo=lo_name, n_hi=len(hi), n_lo=len(lo),
        mean_hi=float(np.mean(hi)), mean_lo=float(np.mean(lo)),
        median_hi=float(np.median(hi)), median_lo=float(np.median(lo)),
        statistic=float(s1), p_one_tailed=float(p1), p_two_tailed=float(p2),
        effect_size_name="matched-pairs rank-biserial r", effect_size=float(r),
        diff_of_means=float(np.mean(d)), ci95_lo=ci_lo, ci95_hi=ci_hi,
        p_two_tailed_floor=p_floor_signrank(len(hi)),
        hedges_g=hedges_g(hi, lo, paired=True),
        g_ci95_lo=g_lo, g_ci95_hi=g_hi,
    ))


def add_permutation(label, family, panel, a, b, a_name, b_name, direction="less"):
    """
    Exact permutation test on the difference in means, enumerating every way of
    splitting the pooled values into groups of the observed sizes. This is what
    the original MP4/MP5 analysis did for n = 4 vs n = 15.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    pooled = np.concatenate([a, b])
    na = len(a)
    obs = a.mean() - b.mean()
    diffs = []
    for combo in combinations(range(len(pooled)), na):
        m = np.zeros(len(pooled), bool)
        m[list(combo)] = True
        diffs.append(pooled[m].mean() - pooled[~m].mean())
    diffs = np.asarray(diffs)
    p_one = float(np.mean(diffs <= obs)) if direction == "less" else float(np.mean(diffs >= obs))
    p_two = float(np.mean(np.abs(diffs) >= abs(obs)))
    ci_lo, ci_hi = boot_ci(a, b)
    g_lo, g_hi = boot_g_ci(a, b)
    rows.append(dict(
        analysis=label, family=family, panel=panel,
        test=f"Exact permutation ({len(diffs)} permutations)",
        group_hi=a_name, group_lo=b_name, n_hi=na, n_lo=len(b),
        mean_hi=float(a.mean()), mean_lo=float(b.mean()),
        median_hi=float(np.median(a)), median_lo=float(np.median(b)),
        statistic=float(obs), p_one_tailed=p_one, p_two_tailed=p_two,
        effect_size_name="difference in means", effect_size=float(obs),
        diff_of_means=float(obs), ci95_lo=ci_lo, ci95_hi=ci_hi,
        p_two_tailed_floor=2.0 / len(diffs),
        hedges_g=hedges_g(a, b),
        g_ci95_lo=g_lo, g_ci95_hi=g_hi,
    ))


# ------------------------------------------------------------------- helpers
def sample_means_from_h5ad(path, gene, phase, group_col, min_cells=MIN_CELLS,
                           site="Stomach", use_raw=True):
    """Sample-level mean expression of `gene`, exactly as the panel scripts do."""
    ad = sc.read_h5ad(path)
    if site is not None and "Sample site" in ad.obs:
        ad = ad[ad.obs["Sample site"] == site]
    ad = ad[ad.obs["Treatment phase"] == phase]
    ad = ad[ad.obs[group_col].isin(["Responsed", "No-response"])].copy()

    src = ad.raw if (use_raw and ad.raw is not None and gene in ad.raw.var_names) else ad
    idx = list(src.var_names).index(gene)
    x = src.X[:, idx]
    x = x.toarray().flatten() if hasattr(x, "toarray") else np.asarray(x).flatten()

    df = pd.DataFrame({
        "sample": ad.obs["sample"].astype(str).values,
        "group": ad.obs[group_col].astype(str).values,
        "value": x,
    })
    keep = df.groupby("sample").size()
    df = df[df["sample"].isin(keep[keep >= min_cells].index)]
    out = df.groupby(["sample", "group"], observed=True)["value"].mean().reset_index()
    return (out.loc[out["group"] == "No-response", "value"].values,
            out.loc[out["group"] == "Responsed", "value"].values)


# =========================================================== 1. CEACAM scRNA
print("[1/9] CEACAM5/6 in pre-treatment epithelium ...")
for gene in ("CEACAM6", "CEACAM5"):
    nr, r = sample_means_from_h5ad(
        EPITHELIAL_H5AD, gene, "Pre", "stomach_pre_grouping")
    add_mw(f"{gene} expression, pre-treatment epithelium",
           "CEACAM", f"Fig 2 ({'N1' if gene == 'CEACAM6' else 'N2'})",
           nr, r, "NR", "R")

# ================================================ 2. CEACAM external (TIGER)
print("[2/9] CEACAM5/6 in PRJEB25780 (TIGER) deconvolved epithelium ...")
epi = pd.read_csv(TIGER_BAYESPRISM_EPI, sep="\t", index_col=0)
meta = pd.read_csv(TIGER_META, sep="\t")
meta = meta[meta["Treatment"] != "Normal"]
common = [s for s in meta["sample_id"] if s in epi.index]
meta = meta[meta["sample_id"].isin(common)].set_index("sample_id")
epi = epi.loc[common]
for gene in ("CEACAM6", "CEACAM5"):
    r = np.log2(epi.loc[meta["response_NR"] == "R", gene].values + 1)
    nr = np.log2(epi.loc[meta["response_NR"] == "N", gene].values + 1)
    add_mw(f"{gene} expression, PRJEB25780 deconvolved epithelium",
           "CEACAM", f"Fig 2 ({'O1' if gene == 'CEACAM6' else 'O2'})",
           nr, r, "NR", "R")

# ============================================================ 3. IHC protein
print("[3/9] CEACAM5/6 immunohistochemistry ...")
ihc = pd.read_csv(PREPARATION / "IHC" / "ceacam_ihc_color_deconv_results.csv")
piv = ihc.pivot_table(index=["patient", "group"], columns="marker",
                      values="staining_pct").reset_index()
piv["combined"] = piv["CEACAM5"] + piv["CEACAM6"]
for col, name in (("CEACAM5", "CEACAM5 only"),
                  ("CEACAM6", "CEACAM6 only"),
                  ("combined", "CEACAM5 + CEACAM6 summed")):
    nr = piv.loc[piv["group"] == "NR", col].values
    r = piv.loc[piv["group"] == "R", col].values
    add_mw(f"IHC staining, {name}", "CEACAM", "Fig 2 (Q)", nr, r, "NR", "R")

# ======================================================= 4. NMF metaprograms
print("[4/9] NMF metaprograms MP4 / MP5 ...")
mp = json.loads((PREPARATION / "Metaprogram_Permutation" / "mp4_permutation_results.json")
                .read_text())
for prog in ("S-MP4", "S-MP5"):
    d = mp[prog]
    pre_r = np.asarray(d["Pre-R_values"], float)
    others = np.concatenate([np.asarray(d[f"{g}_values"], float)
                             for g in ("Post-R", "Pre-NR", "Post-NR")])
    add_permutation(f"{prog} score, pre-treatment R vs all other groups",
                    "Metaprogram", "Fig 2 (K)", pre_r, others,
                    "Pre-R", "All others", direction="less")

# ============================================================== 5. Spatial
print("[5/9] Spatial CEACAM-high vs CEACAM-low regions ...")
spot = pd.read_csv(SPATIAL_SPOT_DATA)
sample_col = "Sample" if "Sample" in spot.columns else "sample"
agg = (spot[spot["CEACAM_group"].isin(["CEACAM-high", "CEACAM-low"])]
       .groupby([sample_col, "CEACAM_group"], observed=True)
       [["neighborhood_epi_density", "distance_to_stroma", "distance_to_immune"]]
       .mean().reset_index())
wide = agg.pivot(index=sample_col, columns="CEACAM_group").dropna()

# The two distance columns are raw Euclidean distances in the x and y of
# spot_data.csv - full-resolution image pixels - and nothing anywhere in that
# path converts them. Until 2026-09-10 this table reported them in those units
# while the Results called the same two numbers micrometres. Author's ruling
# that day: convert both sides, using the project's one derivation of the
# factor, 00_Config/spatial_scale.py (RULES.md rule 5). Fig 3 (H) is a
# proportion and is not touched.
#
# One cohort factor, so this is a pure linear rescale and not a re-weighting:
# the Wilcoxon statistic, both P values, the matched-pairs rank-biserial r,
# Hedges' g and its bootstrap interval are all scale-free and cannot move, and
# the percentile bootstrap interval of the difference of means is the old
# interval times the factor, because boot_ci() resamples INDEX arrays drawn
# from the shared RNG stream and those do not depend on the values. Every other
# row of this table is therefore byte-identical across the change.
# 10_Reproduction/verify_distance_units.py is the check, and it carries the
# mutation that makes it fail.
UM_PER_UNIT = cohort_um_per_unit(spot)
print(f"      distances -> micrometres at {UM_PER_UNIT:.9f} um per array unit")
for _col in ("distance_to_stroma", "distance_to_immune"):
    for _grp in ("CEACAM-high", "CEACAM-low"):
        wide[(_col, _grp)] = wide[(_col, _grp)] * UM_PER_UNIT

for col, label, panel in (
    ("neighborhood_epi_density", "Neighbourhood epithelial density", "Fig 3 (H)"),
    ("distance_to_stroma", "Distance to stroma (µm)", "Fig 3 (I)"),
    ("distance_to_immune", "Distance to immune-rich regions (µm)", "Fig 3 (J)"),
):
    add_wilcoxon(f"{label}, CEACAM-high vs CEACAM-low spots",
                 "Spatial", panel,
                 wide[(col, "CEACAM-high")].values,
                 wide[(col, "CEACAM-low")].values,
                 "CEACAM-high", "CEACAM-low")

# ============================================================== 6. PD-L1
print("[6/9] CD274 (PD-L1) in post-treatment samples ...")
for path, name, panel in (
    (MOMAC_H5AD, "monocytes/macrophages", "Fig 5 (J)"),
    (EPITHELIAL_H5AD, "epithelial cells", "Fig 5 (K)"),
    (FIBROBLAST_H5AD, "fibroblasts", "Fig 5 (L)"),
    (DC_CELLS_H5AD, "dendritic cells", "Fig 5 (M)"),
):
    nr, r = sample_means_from_h5ad(path, "CD274", "Post", "stomach_post_grouping")
    add_mw(f"CD274 expression, post-treatment {name}", "PD-L1", panel,
           nr, r, "NR", "R")


# Hallmark IL-6/JAK/STAT3 signalling, as listed in 05_Figure_5/05_IL6_CD4.
IL6_JAK_STAT3 = [
    "INHBE", "IL17RA", "IRF9", "IL17RB", "MAP3K8", "CCR1", "FAS", "CXCL3", "A2M",
    "CD38", "SOCS3", "TYK2", "GRB2", "CXCL13", "TNFRSF1B", "CXCL1", "CBL", "PF4",
    "CSF1", "IFNGR1", "HMOX1", "TNF", "HAX1", "IL12RB1", "CSF2", "IL2RG", "JUN",
    "ITGA4", "IL18R1", "IL6", "MYD88", "CXCL11", "LEPR", "LTB", "PDGFC", "PTPN11",
    "IFNAR1", "DNTT", "IL1B", "SOCS1", "TNFRSF12A", "PIK3R5", "IL2RA", "CSF2RA",
    "STAT3", "IL13RA1", "BAK1", "TLR2", "CRLF2", "CXCL9", "PIM1", "TNFRSF21",
    "PTPN2", "OSMR", "CSF3R", "IL4R", "IL6ST", "STAM2", "CSF2RB", "EBI3", "STAT2",
    "TNFRSF1A", "IL1R2", "STAT1", "CCL7", "CD14", "TGFB1", "IRF1", "IL3RA",
    "IL10RB", "IL1R1", "CD44", "ITGB3", "ACVRL1", "CXCL10", "IL15RA", "CNTFR",
    "PLA2G2A", "ACVR1B", "IL9R", "LTBR", "CD9", "IFNGR2", "PTPN1", "CD36",
    "REG1A", "IL7",
]

# The pySCENIC regulons plotted in Figure 5D and 5E, as listed in 05_TF.
REGULONS = {
    "BACH1": ["ACSL1", "ARHGAP26", "ASAP1", "AZIN1", "BACH1", "BTG3", "CDC42EP3",
              "CREM", "CSGALNACT2", "DSE", "ELOVL5", "FAM102B", "FAM210A",
              "FNDC3A", "FNDC3B", "GPAT4", "GPCPD1", "HIVEP2", "ITGB1",
              "IVNS1ABP", "JAK1", "JARID2", "KAT6A", "KCNA3", "KDM3A", "KMT2C",
              "KMT2E", "NAB1", "NAMPT", "NFKB1", "NR3C1", "PAG1", "PCNX1",
              "PHF21A", "PIM1", "PPP3CA", "RAP1B", "RASA2", "SEC24A", "SLC44A1",
              "SPEN", "TP53BP2", "TRIP12", "USP12", "ZFYVE16", "ZNF395"],
    "NFKB1": ["ABCA1", "ACSL1", "ACSL4", "AFF4", "AFTPH", "AKT3", "ANKRD12",
              "ARAP2", "ASAP1", "ATP1B3", "ATXN1", "AZIN1", "B3GNT5", "B4GALT5",
              "BACH1", "BASP1", "BAZ1A", "BTG3", "CCNI", "CDC42EP3", "CEP170",
              "CSGALNACT2", "CTNNB1", "CYLD", "DNAJB6", "DSE", "DUSP16", "ELOVL7",
              "EML4", "EPB41L3", "F3", "FAM102B", "FAM107B", "FNBP1", "FNDC3A",
              "FNDC3B", "FRMD6", "GALNTL6", "GPBP1", "HIVEP1", "HIVEP2",
              "IVNS1ABP", "JAK1", "JARID2", "KCNJ2", "KDM7A", "KMT2E", "KPNA4",
              "LDLRAD4", "LYN", "MAPK6", "MIR155HG", "MIR3945HG", "N4BP2",
              "NABP1", "NAMPT", "NFAT5", "NFE2L2", "NFKB1", "NR3C1", "PDE4B",
              "PELI1", "PIM1", "RAP1B", "REL", "SNX9", "SRSF12", "STK26",
              "SUSD6", "TET2", "TP53BP2", "UAP1", "USP12", "WTAP", "ZFYVE16",
              "ZSWIM6"],
}


def sample_scores(adata, genes, min_cells=MIN_CELLS, use_raw=True, ctrl_size=None):
    """
    Sample-level mean module score, exactly as the panel scripts compute it.

    ctrl_size has to match the panel: score_genes picks its control set from
    expression-matched bins, so a different control size gives a different score
    and therefore a different P value. The cytokine panel caps it at 50, the
    regulon panels use the full gene list.
    """
    # Which matrix is the expression is asked once, of shared/expression.py.
    #
    # The line below used to read `adata.raw.var_names if (use_raw and
    # adata.raw is not None) else adata.var_names`, which is guarded - but the
    # score_genes call underneath it was not: it took `use_raw=use_raw`, which
    # defaults to True, and scanpy then dereferences `adata.raw.var_names`
    # itself (scanpy/tools/_score_genes.py:196). Twelve of the thirteen
    # deposited objects carry no .raw, so a reviewer running the capsule
    # against the Zenodo record got an AttributeError here instead of
    # Supplementary Table 6's source table.
    #
    # expression_source() returns .raw where there is one and .X only where the
    # object proves it is the promoted deposit; it raises otherwise, so this is
    # not a fallback. `use_raw` for scanpy is then resolved from that same
    # answer rather than assumed, which is what keeps the two in step: with a
    # .raw present both read .raw, and the numbers do not move.
    src = expression_source(adata) if use_raw else adata
    present = [g for g in genes if g in set(src.var_names)]
    sc.tl.score_genes(adata, gene_list=present, score_name="value",
                      ctrl_size=ctrl_size or min(50, len(present)),
                      use_raw=use_raw and adata.raw is not None)
    df = pd.DataFrame({"sample": adata.obs["sample"].astype(str).values,
                       "group": adata.obs["group"].astype(str).values,
                       "value": adata.obs["value"].values})
    keep = df.groupby("sample").size()
    df = df[df["sample"].isin(keep[keep >= min_cells].index)]
    return df.groupby(["sample", "group"], observed=True)["value"].mean().reset_index()


# ======================================= 7. IL-6/JAK/STAT3 in CD4+ T cells
# Figure 5L. Converted to two-sided with the rest but missing from the first
# version of this sweep, so Table S6 did not list it.
print("[7/9] IL-6/JAK/STAT3 module score in post-treatment CD4+ T cells ...")
ad = sc.read_h5ad(TCD4_H5AD)
ad = ad[ad.obs["Sample site"] == "Stomach"]
ad = ad[ad.obs["Treatment phase"] == "Post"].copy()
ad.obs["group"] = ad.obs["stomach_post_grouping"].astype(str).map(
    {"Responsed": "R", "No-response": "NR"}).astype(str)
ad = ad[ad.obs["group"].isin(["R", "NR"])].copy()
sm = sample_scores(ad, IL6_JAK_STAT3)
add_mw("IL-6/JAK/STAT3 module score, post-treatment CD4+ T cells",
       "Cytokine", "Fig 5 (IL6-CD4)",
       sm.loc[sm["group"] == "NR", "value"].values,
       sm.loc[sm["group"] == "R", "value"].values, "NR", "R")

# ==================================== 8. NF-kB regulon activity in IL-1b+ Mac
# Figures 5D and 5E, same omission.
print("[8/9] BACH1 and NFKB1 regulon activity in IL-1beta+ macrophages ...")
ad = sc.read_h5ad(MOMAC_H5AD)
if ad.raw is not None:
    ad = ad.raw.to_adata()
ad = ad[ad.obs["minor_cell_state"].astype(str).str.contains("C3")].copy()
pre = ad.obs["stomach_pre_grouping"].astype(str)
post = ad.obs["stomach_post_grouping"].astype(str)
ad.obs["group"] = "Other"
ad.obs.loc[pre == "Responsed", "group"] = "Pre-R"
ad.obs.loc[pre == "No-response", "group"] = "Pre-NR"
ad.obs.loc[post == "Responsed", "group"] = "Post-R"
ad.obs.loc[post == "No-response", "group"] = "Post-NR"
# Scored over every C3 cell, including the ones that fall outside the four
# groups, because that is the cell set the panel scores on and score_genes picks
# its control genes from expression bins computed over whatever is present.
for tf, genes in REGULONS.items():
    present = [g for g in genes if g in set(ad.var_names)]
    sm = sample_scores(ad, genes, min_cells=1, use_raw=False,
                       ctrl_size=len(present))
    sm = sm[sm["group"] != "Other"]
    post_r = sm.loc[sm["group"] == "Post-R", "value"].values
    others = sm.loc[sm["group"] != "Post-R", "value"].values
    add_permutation(f"{tf} regulon activity, post-treatment R vs all other groups",
                    "Regulon", f"Fig 5 (regulon-{tf})",
                    post_r, others, "Post-R", "All others", direction="less")


# ================================ 9. CEACAM5/6+ proportion, tumor-content adjusted
# Supplementary Figure S2D. The submitted panel printed the one-sided exact
# value (P = 0.03) and was missed by the first version of this sweep, so Table
# S6 did not list it. Computed exactly as
# 03_Final_Panels/02_Figure_2/02_E2/create_c2_proportion_adjusted.py
# does: the per-sample percentage of C2_Epi_CEACAM6 cells among pre-treatment
# gastric epithelial cells of the response-classified samples, regressed by
# ordinary least squares on the sample's mean tumor score, and the residuals
# compared NR vs R. It comes last so that the shared bootstrap stream seeded
# above reaches every earlier row exactly as it did before this one existed.
# The tumour score is an obs column of the epithelial object, beside the
# cell-state and grouping columns the proportion is taken over, so one read
# supplies the whole frame.
print("[9/9] CEACAM5/6+ proportion adjusted for tumor content (Fig. S2D) ...")
obs = sc.read_h5ad(EPITHELIAL_DEPOSIT_H5AD, backed="r").obs.copy()
keep = ((obs["Treatment phase"] == "Pre")
        & obs["stomach_pre_grouping"].isin(["Responsed", "No-response"]))
sub = obs[keep].copy()
sub["is_C2"] = (sub["minor_cell_state"] == "C2_Epi_CEACAM6").astype(int)
per = sub.groupby("sample", observed=True).agg(
    C2_proportion=("is_C2", lambda x: x.mean() * 100),
    mean_tumor_score=("tumor_score", "mean"),
    group=("stomach_pre_grouping", "first"),
).reset_index().dropna()
slope, intercept, _, _, _ = stats.linregress(per["mean_tumor_score"].values,
                                             per["C2_proportion"].values)
per["residual"] = (per["C2_proportion"].values
                   - (slope * per["mean_tumor_score"].values + intercept))
add_mw("CEACAM5/6+ proportion, tumor-content adjusted (Fig. S2D)",
       "CEACAM", "Fig S2 (D)",
       per.loc[per["group"] == "No-response", "residual"].values,
       per.loc[per["group"] == "Responsed", "residual"].values, "NR", "R")


# ---------------------------------------------------------------- assemble
df = pd.DataFrame(rows)

# Multiple-testing correction within each hypothesis family, on the two-sided P.
df["p_two_tailed_BH"] = np.nan
for fam, grp in df.groupby("family"):
    df.loc[grp.index, "p_two_tailed_BH"] = multipletests(
        grp["p_two_tailed"], method="fdr_bh")[1]

df["sig_one_tailed"] = df["p_one_tailed"] < 0.05
df["sig_two_tailed"] = df["p_two_tailed"] < 0.05
df["sig_two_tailed_BH"] = df["p_two_tailed_BH"] < 0.05
df["ci_excludes_zero"] = (df["ci95_lo"] > 0) | (df["ci95_hi"] < 0)
# A test is at its resolution limit when the observed two-sided P is within one
# step of the smallest value the exact null can produce at these group sizes.
df["at_resolution_limit"] = df["p_two_tailed"] <= 2.5 * df["p_two_tailed_floor"]


def verdict(r):
    if r["sig_two_tailed"]:
        return "Robust: significant two-sided"
    if r["sig_one_tailed"] and r["ci_excludes_zero"]:
        return "Weakened: not significant two-sided, but CI excludes zero"
    if r["sig_one_tailed"]:
        return "Lost: significant one-tailed only"
    return "Not significant either way"


df["verdict"] = df.apply(verdict, axis=1)
df = df.sort_values(["family", "panel", "analysis"])
df.to_csv(OUT / "twosided_sweep.csv", index=False)

# ------------------------------------------------------------------ report
lines = ["TWO-SIDED STATISTICS SWEEP - Reviewer 1 point R1.3c", "=" * 100, ""]
for fam, grp in df.groupby("family"):
    lines.append(f"[{fam}]")
    for _, r in grp.iterrows():
        lines.append(f"  {r['analysis']}   ({r['panel']}, {r['test']})")
        lines.append(
            f"     {r['group_hi']} n={r['n_hi']} mean={r['mean_hi']:.4g}   "
            f"{r['group_lo']} n={r['n_lo']} mean={r['mean_lo']:.4g}")
        lines.append(
            f"     P one-tailed = {r['p_one_tailed']:.4g}   "
            f"P two-tailed = {r['p_two_tailed']:.4g}   "
            f"BH within family = {r['p_two_tailed_BH']:.4g}")
        lines.append(
            f"     smallest two-sided P attainable at n={r['n_hi']} vs "
            f"{r['n_lo']} = {r['p_two_tailed_floor']:.4g}"
            + ("   <-- AT THE RESOLUTION LIMIT OF THE TEST"
               if r["at_resolution_limit"] else ""))
        lines.append(
            f"     {r['effect_size_name']} = {r['effect_size']:.3f}   "
            f"diff of means = {r['diff_of_means']:.4g} "
            f"[95% CI {r['ci95_lo']:.4g}, {r['ci95_hi']:.4g}]")
        lines.append(f"     -> {r['verdict']}")
        lines.append("")
    lines.append("")

lines.append("SUMMARY")
lines.append("-" * 100)
lines.append(f"  Comparisons rerun                              : {len(df)}")
lines.append(f"  Significant one-tailed (as published)          : {int(df['sig_one_tailed'].sum())}")
lines.append(f"  Still significant two-sided                    : {int(df['sig_two_tailed'].sum())}")
lines.append(f"  Still significant after BH within family       : {int(df['sig_two_tailed_BH'].sum())}")
lines.append(f"  Bootstrap 95% CI excludes zero                 : {int(df['ci_excludes_zero'].sum())}")
lines.append(f"  At the resolution limit of the exact null      : {int(df['at_resolution_limit'].sum())}")
lines.append("")
for v, n in df["verdict"].value_counts().items():
    lines.append(f"    {v}: {n}")

report = "\n".join(lines)
(OUT / "twosided_sweep_report.txt").write_text(report, encoding="utf-8")
print("\n" + report)
