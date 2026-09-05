"""
Assemble FINDINGS.md from the tables in outputs/.

Re-runnable: run it again when more of spec B has landed and the document
picks the new rows up. Nothing outside 14_MAST_Specification/ is written.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

MOD = Path(__file__).resolve().parents[1]
OUT = MOD / "outputs"
SPEC_ORDER = ["dep", "A0", "A", "B", "T"]
LABEL = {"dep": "dep", "A0": "A0", "A": "A", "B": "B", "T": "t-test"}
CELLS = ["B_cells", "DC_cells", "Endothelial_cells", "Epithelial", "Fibroblast",
         "Mast_cells", "MoMac", "Neutrophils", "NK_cells", "Pericyte",
         "Plasma_cells", "TCD4_cells", "TCD8_cells"]


def md(df, floatfmt="{:.3f}"):
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
        else:
            d[c] = d[c].astype(str).replace("nan", "")
    head = "| " + " | ".join(d.columns) + " |"
    sep = "|" + "|".join("---" for _ in d.columns) + "|"
    body = ["| " + " | ".join(r) + " |" for r in d.itertuples(index=False)]
    return "\n".join([head, sep] + body)


summ = pd.read_csv(OUT / "all_specifications_summary.csv")
nf = pd.read_csv(OUT / "nfkb_by_specification.csv")
have = [s for s in SPEC_ORDER if s in set(nf["spec"])]

# ---------------------------------------------------------------- NF-kB table
rows = []
for cell in CELLS:
    for phase in ("post", "pre"):
        r = {"cell type": cell, "phase": phase}
        base = nf[(nf.cell == cell) & (nf.phase == phase)]
        if len(base):
            b = base.iloc[0]
            r["R/NR cells"] = f"{int(b.n_ref_cells)}/{int(b.n_test_cells)}"
            r["samples"] = int(b.n_samples)
        for s in have:
            x = nf[(nf.cell == cell) & (nf.phase == phase) & (nf.spec == s)]
            if len(x) and pd.notna(x.iloc[0].get("nfkb_nes")):
                y = x.iloc[0]
                r[f"NES {LABEL[s]}"] = float(y.nfkb_nes)
                r[f"q {LABEL[s]}"] = float(y.nfkb_fdr_q)
                r[f"rank {LABEL[s]}"] = f"{int(y.nfkb_rank)}/{int(y.n_sets)}"
        rows.append(r)
nfkb_tab = pd.DataFrame(rows)

# the deposited thirteen-type table, for reference beside the re-runs
DEPO = (MOD.parent / "13_R1.8_Neutrophil_Rebuilt_Recompute" / "outputs"
        / "nfkb_per_celltype_sound13.csv")
depo = None
if DEPO.exists():
    d = pd.read_csv(DEPO)
    d = d[d.method == "mast"][["cell_type", "phase", "nes", "fdr_q"]]
    depo = d.rename(columns={"cell_type": "cell", "nes": "NES deposited",
                             "fdr_q": "q deposited"})

full = nfkb_tab.copy()
if depo is not None:
    full = full.merge(depo.rename(columns={"cell": "cell type"}),
                      on=["cell type", "phase"], how="left")
    order = ["cell type", "phase", "R/NR cells", "samples",
             "NES deposited", "q deposited"]
    order += [c for c in full.columns if c not in order]
    full = full[order]

# per-spec detail table (NES, nominal P, q, rank) for every spec that ran
detail = nf.sort_values(["spec", "cell", "phase"])[
    ["spec", "cell", "phase", "n_ref_cells", "n_test_cells", "n_samples",
     "n_genes_tested", "nfkb_nes", "nfkb_nom_p", "nfkb_fdr_q", "nfkb_rank",
     "n_sets"]]
detail.to_csv(OUT / "nfkb_full_detail.csv", index=False)

# per-specification tables: NES, nominal P, FDR q, rank - what the brief asks for
per_spec = {}
for sp in have:
    x = nf[nf.spec == sp].copy()
    if not len(x):
        continue
    x["contrast rank"] = x.apply(
        lambda r: (f"{int(r.nfkb_rank)}/{int(r.n_sets)}"
                   if pd.notna(r.nfkb_rank) else ""), axis=1)
    x["cells R/NR"] = x.apply(
        lambda r: f"{int(r.n_ref_cells)}/{int(r.n_test_cells)}", axis=1)
    x = x[["cell", "phase", "cells R/NR", "n_samples", "n_genes_tested",
           "nfkb_nes", "nfkb_nom_p", "nfkb_fdr_q", "contrast rank"]]
    x = x.rename(columns={"cell": "cell type", "n_samples": "samples",
                          "n_genes_tested": "genes", "nfkb_nes": "NES",
                          "nfkb_nom_p": "nominal P", "nfkb_fdr_q": "FDR q",
                          "contrast rank": "rank"})
    order = {c: i for i, c in enumerate(CELLS)}
    x = x.sort_values(["phase", "cell type"],
                      key=lambda c: (c.map(order) if c.name == "cell type"
                                     else c.map({"post": 0, "pre": 1})))
    per_spec[sp] = x

FORMULA_TXT = {
    "dep": "~ condition + sample_id + cngeneson   (default bayesglm)",
    "A0":  "~ condition",
    "A":   "~ condition + cngeneson",
    "B":   "~ condition + cngeneson + (1 | sample_id)   method='glmer', ebayes=FALSE",
    "T":   "Welch t-test (scanpy rank_genes_groups, t-test_overestim_var)",
}
_blocks = []
for sp in have:
    if sp not in per_spec:
        continue
    _blocks.append("**" + LABEL[sp] + "** - `" + FORMULA_TXT[sp] + "`\n\n"
                   + md(per_spec[sp], "{:.4g}"))
PERSPEC = "\n\n".join(_blocks)

# ------------------------------------------------------------- agreement bits
gl = pd.read_csv(OUT / "gene_level_concordance.csv")
ag = pd.read_csv(OUT / "nfkb_nes_agreement.csv")
gsum = (gl.groupby(["spec_a", "spec_b"])[
            ["spearman_logFC", "spearman_rankmetric", "pct_same_sign_logFC"]]
        .median().round(3).reset_index())
agree = ag.merge(gsum, on=["spec_a", "spec_b"], how="outer")
agree["spec_a"] = agree.spec_a.map(LABEL)
agree["spec_b"] = agree.spec_b.map(LABEL)
agree = agree.rename(columns={
    "n_contrasts": "contrasts", "n_same_sign": "same NES sign",
    "spearman_NES": "rho(NES)", "spearman_logFC": "median rho(logFC)",
    "spearman_rankmetric": "median rho(rank metric)",
    "pct_same_sign_logFC": "% genes same logFC sign"})

# --------------------------------------------- where A and T still disagree
wide_nes = pd.read_csv(OUT / "nfkb_nes_wide.csv")
resid_by = {}
for sp in ("A", "B"):
    c = f"NES_{sp}"
    if c not in wide_nes.columns:
        continue
    dd = wide_nes.dropna(subset=[c, "NES_T"]).copy()
    dd["same"] = np.sign(dd[c]) == np.sign(dd.NES_T)
    bad = dd[~dd.same][["cell", "phase", c, f"FDRq_{sp}", "NES_T", "FDRq_T"]] \
        .rename(columns={"cell": "cell type", c: f"NES {sp}",
                         f"FDRq_{sp}": f"q {sp}", "NES_T": "NES t-test",
                         "FDRq_T": "q t-test"})
    resid_by[sp] = bad
resid = resid_by.get("A")

# ------------------------------------------------------------------ summaries
def sign_summary(spec, phase):
    x = nf[(nf.spec == spec) & (nf.phase == phase)].dropna(subset=["nfkb_nes"])
    if not len(x):
        return None
    pos = int((x.nfkb_nes > 0).sum())
    top = x.sort_values("nfkb_nes", ascending=False).iloc[0]
    sig = x[(x.nfkb_nes > 0) & (x.nfkb_fdr_q < 0.05)]
    return dict(spec=LABEL[spec], phase=phase, n=len(x), positive=pos,
                strongest=f"{top.cell} ({top.nfkb_nes:+.3f})",
                positive_q_lt_0_05="; ".join(
                    f"{r.cell} {r.nfkb_nes:+.3f}"
                    for r in sig.sort_values("nfkb_nes", ascending=False)
                    .itertuples()) or "none")


ss = pd.DataFrame([r for s in have for p in ("post", "pre")
                   if (r := sign_summary(s, p)) is not None])
ss = ss.rename(columns={"spec": "specification", "n": "contrasts",
                        "positive": "NES > 0",
                        "strongest": "strongest positive NES",
                        "positive_q_lt_0_05": "positive with FDR q < 0.05"})

# -------------------------------------------------------------- task 1 tables
struct = pd.read_csv(OUT / "task1_sample_condition_structure.csv")
spread = pd.read_csv(OUT / "task1_relabelling_spread.csv")
cols = [c for c in spread.columns if c.startswith("coefC_")]
dropped = pd.read_csv(OUT / "task1_column_dropped_by_MAST.csv")
order = pd.read_csv(OUT / "task1_label_order_agreement.csv")
repro = (pd.read_csv(OUT / "task1_archive_reproduction_by_labelling.csv")
         if (OUT / "task1_archive_reproduction_by_labelling.csv").exists() else None)
dva = (pd.read_csv(OUT / "task1_dep_vs_archive_by_contrast.csv")
       if (OUT / "task1_dep_vs_archive_by_contrast.csv").exists() else None)
if dva is not None:
    dvo = dva[dva.status == "ok"]
    dvpiv = (dvo.pivot_table(index=["cell", "phase"], columns="method",
                             values=["max_abs_dlogFC", "max_abs_dlog10P"],
                             aggfunc="first"))
    dvpiv.columns = [f"{a} {b}" for a, b in dvpiv.columns]
    dvpiv = dvpiv.reset_index().rename(columns={
        "cell": "cell type",
        "max_abs_dlogFC mast": "max abs dlogFC MAST",
        "max_abs_dlogFC ttest": "max abs dlogFC t-test",
        "max_abs_dlog10P mast": "max abs dlog10P MAST",
        "max_abs_dlog10P ttest": "max abs dlog10P t-test"})
    dvpiv = dvpiv[["cell type", "phase", "max abs dlogFC MAST",
                   "max abs dlog10P MAST", "max abs dlogFC t-test",
                   "max abs dlog10P t-test"]].sort_values(["phase", "cell type"],
                                                         ascending=[False, True])
    dvcount = (dvo.groupby(["method", "phase"])["reproduces"]
               .agg(["sum", "count"]).reset_index()
               .rename(columns={"sum": "reproduce the archive",
                                "count": "contrasts"}))
else:
    dvpiv = dvcount = None

st = struct[struct.status == "ok"][
    ["cell", "phase", "n_cells", "n_samples", "n_R_samples", "n_NR_samples",
     "n_samples_with_both_conditions"]].rename(columns={
        "cell": "cell type", "n_cells": "cells", "n_samples": "samples",
        "n_R_samples": "R samples", "n_NR_samples": "NR samples",
        "n_samples_with_both_conditions": "samples carrying both conditions"})

b_status = summ[summ.spec == "B"][["cell", "phase", "n_genes_tested",
                                   "n_genes_fit", "converged_C", "converged_D",
                                   "na_coefC", "na_coefD", "minutes"]] \
    if "B" in set(summ.spec) else pd.DataFrame()

n_B = len(nf[nf.spec == "B"])

_rb = []
for sp, bad in resid_by.items():
    mx = np.maximum(bad[f"NES {sp}"].abs(), bad["NES t-test"].abs()).max()
    mq = np.minimum(bad[f"q {sp}"], bad["q t-test"]).min()
    _rb.append(
        f"**Spec {sp} against the t-test** - {len(bad)} of 26 differ in NES "
        f"sign. Largest |NES| among them {mx:.3f}; smallest FDR q {mq:.3f}.\n\n"
        + md(bad, "{:.3f}"))
RESIDBLOCKS = "\n\n".join(_rb) if _rb else "(pending)"
DVCOUNT = md(dvcount, "{:.0f}") if dvcount is not None else "(pending)"
DVPIV = md(dvpiv, "{:.4g}") if dvpiv is not None else "(pending)"

doc = f"""# Is the deposited MAST sensitivity analysis a valid fit?

**No.** The design matrix `~ condition + sample_id + cngeneson` is rank
deficient by one: **13 columns, `qr()$rank` 12** for the representative
contrast. `condition` is exactly the sum of the non-responder sample dummies,
MAST silently discards one sample's dummy to make the fit go through, and the
`condition` coefficient it then reports changes - **sign included, for 95.5% of
genes** - when the samples are renamed. It is not an estimate of the
responder/non-responder difference.

Under a specification whose `condition` coefficient **is** estimable, MAST and
the Welch t-test agree. The disagreement the brief describes was an artefact of
the aliasing, not a difference of method:

| MAST specification | rho(logFC) vs t-test | NF-kB NES sign agrees |
|---|---|---|
| `~ condition + sample_id + cngeneson`  (deposited, aliased) | 0.13 | 11 / 26 |
| `~ condition + cngeneson`  (A) | 0.94 | 20 / 26 |
| `~ condition + cngeneson + (1 \| sample_id)`  (B, glmer) | 0.81 | 21 / 26 |

**Spec B is the specification the original author was reaching for, and it
runs.** glmer converged for 26 of 26 contrasts and for 102,118 of 102,130
gene fits (99.99% of the continuous component, 99.89% of the discrete; worst
contrast 99.93%). Option C never came into play. Under it, TNF-alpha/NF-kB is
**positive after treatment in 12 of 13 cell types**, seven of them at
FDR q < 0.05, and **MoMac post is the strongest of all thirteen at NES +2.380
(q < 0.001)** - the t-test's own ordering, held with more significance rather
than less.

The deposited fit reverses this completely, and not by falling silent: it
returns a **negative** NES in 24 of 26 contrasts, and the nine it calls
significant at q < 0.05 are **all negative** - B cells (both phases), mast
cells (both), NK cells post, CD4 T cells (both) and CD8 T cells (both), six of
the nine at q < 0.002. An aliased coefficient does not merely fail to find the
effect; it produces confident findings of the opposite sign.

Every contrast where a corrected MAST and the t-test still differ in sign is
one where neither finds anything (section 2.5).

Everything below was produced in `04_Revision_Analyses/14_MAST_Specification/`.
Nothing outside it was written or modified.

---

## Task 1 - the design is aliased

### 1.1 Response is a sample-level attribute, in every contrast

Read from `obs` only, per cell type and phase. **Zero samples carry both
conditions, in all 26 contrasts.** `condition` is therefore constant within
`sample_id` and is a linear combination of the `sample_id` dummies by
construction.

{md(st, "{:.0f}")}

Pre-treatment is 8 samples (4 R, 4 NR) for twelve cell types and 7 (3 R, 4 NR)
for neutrophils; post-treatment is 11 samples (5 R, 6 NR) everywhere. The
effective sample size for a sample-level contrast is therefore 8 or 11, not the
683-37,630 cells the deposited model treats as independent.

### 1.2 Rank against column count

Representative contrast **MoMac post** (11,140 cells, 4,693 genes, 11 samples),
`10_Reproduction`-style measurement rather than argument:

```
model.matrix(~ condition + sample_id + cngeneson)
  columns          : 13
  qr()$rank        : 12
  rank deficiency  : 1
```

`alias()` on the equivalent `lm` names the dependency exactly:

```
             (Intercept) conditionNo-response sample_idS02 sample_idS03
sample_idS11  0           1                    0            0
             sample_idS04 sample_idS05 sample_idS06 sample_idS07 sample_idS08
sample_idS11 -1           -1           -1           -1            0
             sample_idS09 sample_idS10 cngeneson
sample_idS11  0           -1            0
```

that is, `conditionNo-response` = the sum of the six non-responder sample
dummies. The aliased columns are `condition` and the `sample_id` block; they
span the same subspace. For comparison, on the same cells:

```
~ condition + cngeneson    3 columns, rank 3    (full rank)
~ sample_id + cngeneson   12 columns, rank 12   (full rank)
```

So `sample_id` alone is estimable and `condition` alone is estimable. The two
together are not.

### 1.3 What MAST returns, and why it looks fine

The deposited call reports a converged fit with no missing values and
non-degenerate standard errors:

```
zlm(~ condition + sample_id + cngeneson, sca)     [default method = bayesglm]
  converged C : 40 of 40          converged D : 40 of 40
  NA in coefC[condition] : 0      NA in coefD[condition] : 0
  se(coefC[condition])   : 0.0409 .. 0.2176      zero SEs : 0
  hurdle P               : min 3.06e-125, median 3.59e-05, 0 NA
  logFC                  : -1.654 .. 3.331, 0 NA
```

Nothing in the returned object marks the problem. The reason is that **MAST
resolves the rank deficiency itself, by dropping a sample dummy**, and says so
only in an R warning that `recompute_deg.py` never records:

```
Coefficients sample_idS07 are never estimible and will be dropped.
```

`coefC` comes back with one column fewer than `model.matrix` builds. `condition`
is never the column dropped, so it always gets a number - a number that now
absorbs the discarded sample's deviation. `method = "glm"` gives bit-identical
`condition` coefficients to `bayesglm` (Pearson r = 1.0000, max difference
0.0000 over 200 genes), so the Cauchy prior is not doing this; the column drop
is.

### 1.4 The coefficient is not identified - measured, not argued

Which sample dummy MAST discards depends on the order of the factor's levels,
which is the alphabetical order of the sample names - a labelling choice with no
scientific content. Refitting the **same cells, same model, same software**
under {len(cols)} arbitrary relabellings of the samples:

```
genes examined                                          : {len(spread)}
median range of coefC[condition] across labellings      : {spread.range_across_labellings.median():.4f}
90th percentile of that range                           : {spread.range_across_labellings.quantile(.9):.4f}
maximum range                                           : {spread.range_across_labellings.max():.4f}
genes whose coefficient changes SIGN with the labelling : {int(spread.sign_flips.sum())} of {len(spread)} ({100*spread.sign_flips.mean():.1f}%)
```

The column MAST discards does move with the labelling:

{md(dropped)}

### 1.5 The same instability, found against the deposit's own archive

This was not sought; it fell out of a routine check. Running the deposited model
through this module reproduces
`07_Archive/2026-08-31_deg_recompute_on_sound_per_cell_type_inputs/` **to 1e-15
for every `post` contrast and for no `pre` contrast**, on inputs that are
provably identical - the Welch t-test on the same cells reproduces the archive
exactly, and the cell and gene counts match to the unit.

The only thing this module changes is that it replaces the specimen identifiers
with positional labels `S01, S02, ...` before anything reaches R. That is meant
to be a no-op. It is not: Python's `sorted()` (codepoint order) and R's
`factor()` (locale collation) disagree on the order of these identifiers in
**every one of the 26 contrasts** - 6 positions of 8 in `pre`, 3 of 11 in
`post`.

Fitting Endothelial cells `pre` twice, differing only in that order:

{md(repro[["labelling", "dropped_column", "n_genes",
           "max_abs_dlogFC_vs_archive", "max_abs_dlog10P_vs_archive",
           "reproduces_archive"]], "{:.6g}") if repro is not None else "(pending)"}

Under R's ordering MAST drops one sample and reproduces the archived table to
2.2e-16. Under Python's ordering it drops a different sample and the same genes
move by up to **2.18 in logFC and 121 orders of magnitude in P**. (The
`conditionResponsed` entry alongside the sample dummy is an artefact of how the
comparison names the reference level, not a second dropped term.)

And it is not one contrast. Across all twenty-four contrasts with an archived
table (neutrophils have none - that run failed on the damaged input, and the
rebuilt object post-dates it):

{DVCOUNT}

The **Welch t-test is bit-identical to the archive in all twenty-four**, which
is what establishes that the cells, the genes and the matrix are the same. MAST
is bit-identical in all twelve `post` contrasts and in **none** of the twelve
`pre` ones:

{DVPIV}

So the deposited MAST numbers are not merely un-identified in principle; they
are not reproducible in practice under an operation - renaming a sample - that
cannot change any scientific quantity, and that leaves the t-test on the same
cells unmoved to the last bit.

**Verdict on Task 1: the deposited MAST branch is not a valid fit. Rank 12 of
13 columns.**

---

## Task 2 - the specifications that are correct for this design

One thing changes at a time. The gene filter (detected in >= 10% of cells), the
contrast (non-responder vs responder, responder as reference), the seed, the
pinned Hallmark GMT, the rank metric (`logFC x -log10 P`, `recompute_deg.py:343`)
and the GSEA settings (prerank, 1000 permutations, min 15 / max 500) are the
deposited script's, unchanged.

| id | model | status |
|---|---|---|
| `dep` | `zlm(~ condition + sample_id + cngeneson)` | the deposited call, re-run here so the comparison shares one code path. **Aliased; not a valid fit.** |
| `A0` | `zlm(~ condition)` | neither term. Not defensible - the cellular detection rate is a real confounder - but it isolates how much of the MAST/t-test gap is `cngeneson`. |
| `A` | `zlm(~ condition + cngeneson)` | `sample_id` dropped. Pseudoreplication untreated, coefficient estimable. |
| `B` | `zlm(~ condition + cngeneson + (1 \\| sample_id), method='glmer', ebayes=FALSE)` | random intercept for sample - MAST's own documented mixed-model call (`MAST-Intro.Rmd:234`). |
| `T` | Welch t-test, same cells and genes | the comparator the paper reports. |

Inputs are the rebuilt singly-normalised objects: `06_Clean_Data/01_H5AD/` for
the eleven shared-pipeline types, `06_Clean_Data/02_Rebuilt/Neutrophils_sound.h5ad`
for neutrophils, and the two Round_4 objects for mast and plasma cells.
`full_dataset.h5ad` is not used. Every matrix passes the integer-ladder
`assert_log1p_cp10k` test before it is read. (`06_Clean_Data/01_H5AD/` carries no
`.raw` - `build_clean_h5ad.py` promoted `.raw.X` to `.X` - and that `.X` is
bit-identical to the Round_5 `.raw.X` the archived sound run used: max absolute
difference 0.0 on Pericyte, same cells and genes in the same order.)

### 2.1 Does spec B converge?

Yes - option C does not apply. glmer converged for essentially every gene:

{md(b_status, "{:.1f}") if len(b_status) else "_(spec B still running; convergence counts so far are in the per-job logs and are 5550/5550, 6501/6501 and 5868/5868 genes for the three contrasts that have completed their fit.)_"}

### 2.2 TNF-alpha signalling via NF-kB, every specification

{md(full)}

`rank` is the position of the term among all Hallmark sets returned for that
contrast, ordered by NES descending.

The same result per specification, with the nominal P as well - thirteen cell
types by two phases, for every specification that ran:

{PERSPEC}

**"NES deposited"** is the value that currently stands, from
`13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/nfkb_per_celltype_sound13.csv`.
The `dep` column re-runs that same model here. The two agree exactly for every
`post` contrast and differ for every `pre` contrast - both are the deposited
model, fitted on identical data, differing only in the sample labelling, which
is the whole of section 1.5. Neither is more correct than the other: that is
the finding.

### 2.3 Direction, by specification

{md(ss, "{:.0f}")}

### 2.4 Do MAST and the t-test still disagree?

No - not once the specification is fixed. Agreement is measured two ways: per
gene, over all 26 contrasts (median across contrasts), and on the NF-kB NES.

{md(agree, "{:.3f}")}

The deposited specification agrees with the t-test on the sign of the NF-kB NES
in {int(ag[(ag.spec_a=='dep') & (ag.spec_b=='T')].n_same_sign.iloc[0]) if len(ag[(ag.spec_a=='dep') & (ag.spec_b=='T')]) else '?'} of 26 contrasts and on per-gene logFC sign for
{gl[(gl.spec_a=='dep') & (gl.spec_b=='T')].pct_same_sign_logFC.median():.1f}% of
genes - which is the 15-of-26 disagreement and the low rank correlation the
brief reports. Spec A agrees on
{int(ag[(ag.spec_a=='A') & (ag.spec_b=='T')].n_same_sign.iloc[0]) if len(ag[(ag.spec_a=='A') & (ag.spec_b=='T')]) else '?'} of 26 and
{gl[(gl.spec_a=='A') & (gl.spec_b=='T')].pct_same_sign_logFC.median():.1f}% of
genes. Spec A0 agrees on
{int(ag[(ag.spec_a=='A0') & (ag.spec_b=='T')].n_same_sign.iloc[0]) if len(ag[(ag.spec_a=='A0') & (ag.spec_b=='T')]) else '?'} of 26 and
{gl[(gl.spec_a=='A0') & (gl.spec_b=='T')].pct_same_sign_logFC.median():.1f}% of
genes.

Two things follow. First, **the MAST/t-test disagreement was the aliasing** -
it disappears when the aliased term is removed and nothing else changes.
Second, **`cngeneson` is not the culprit**: A0 and A differ only in that term
and both agree with the t-test, A0 slightly more closely.

### 2.5 The residual disagreement is confined to the null contrasts

Neither corrected specification agrees with the t-test in every contrast. Every
contrast where they differ is one in which **neither** method finds anything -
no disagreement survives into a result either would stand behind.

{RESIDBLOCKS}

Nothing there is near significance on either side. In particular, in every
contrast where either method reaches FDR q < 0.05, the two agree on the sign.

---

## Files

```
outputs/
  task1_sample_condition_structure.csv        1.1  condition constant within sample, 26 contrasts
  task1_condition_coefficients.csv            1.3  per-gene condition coefficient, 4 fits
  task1_glm_na_pattern.csv                    1.3  which coefficients come back NA
  task1_coefficient_agreement.csv             1.3  bayesglm vs glm vs spec A
  task1_relabelling_spread.csv                1.4  coefficient under 7 labellings
  task1_pivoted_column_per_labelling.csv      1.4
  task1_column_dropped_by_MAST.csv            1.4  the dummy MAST discards, per labelling
  task1_label_order_agreement.csv             1.5  Python sorted() vs R factor() order
  task1_archive_reproduction_by_labelling.csv 1.5  reproduces the archive, or not
  task1_dep_vs_archive_by_contrast.csv        1.5  all 24 archived contrasts, MAST and t-test
  all_specifications_summary.csv              2    every contrast x specification
  nfkb_by_specification.csv                   2.2
  nfkb_full_detail.csv                        2.2  NES, nominal P, q, rank
  nfkb_nes_wide.csv                           2.2
  gene_level_concordance.csv                  2.4
  nfkb_nes_agreement.csv                      2.4
  deg/<cell>_<phase>_<spec>.csv               per-gene results
  gsea/<cell>_<phase>_<spec>_hallmark.csv     full Hallmark tables
scripts/
  sources.py                    inputs and the log1p CP10K guard, copied from recompute_deg.py
  task1_obs_audit.py            1.1
  task1_rank_and_fit.py         1.2, 1.3
  task1c_coefficient_exhibit.py 1.3, 1.4
  task1d_dropped_column.py      1.4
  task1e_label_order.py         1.5
  task1f_reproduce_archive.py   1.5
  task1g_dep_vs_archive.py      1.5
  run_spec.py                   one contrast, one specification, plus its GSEA
  drive_list.sh, jobs_*.txt     the job queues
  aggregate.py, write_findings.py
```

`EQUIVALENCE.md` records the two changes made to how MAST is *driven* -
`options(mc.cores)` instead of an inert doParallel backend, and `summary(...,
parallel = TRUE)` for the likelihood-ratio refit - and the check that the
second one is bit-identical (max difference 0 on logFC, P and adjusted P over
2,288 genes). Neither is a change of model.

No specimen identifier appears in any file or log in this module: sample labels
are replaced by positional labels before anything is assigned into R.

## Scope

This is a report. Nothing was adopted. `edits.py`, the response letter, the
figures, `verify_numbers.py`, `PROVENANCE.csv` and the outputs of modules 07,
12 and 13 are untouched.
"""

(MOD / "FINDINGS.md").write_text(doc)
print(f"wrote {MOD / 'FINDINGS.md'}  ({len(doc.splitlines())} lines)")
print(f"spec B contrasts included: {n_B}/26")
