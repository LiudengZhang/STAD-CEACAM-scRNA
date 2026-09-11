> **Upstream record.** Paths in this document refer to the upstream
> processing pipeline that produced the deposited intermediates; it is not
> part of this release's run path. The files under `upstream/` are included
> as a record of how the inputs were produced and are not executed by
> `./run` or `_run_all_panels.sh`.

# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# BayesPrism_TCGA - where this code came from

Source: `submission-tree/temp-workspace/02162026_BayesPrism`

Produces: epithelial expression deconvolved from TCGA-STAD bulk

Consumed by: printed Figure 3B and Figure 3E (build directories `03_C` and
`03_D` - see the correction at the end of this file), the TCGA survival
analyses, and R1.6 immune exclusion.

| file | md5 | bytes |
|---|---|---|
| `step1_prepare_tcga_bulk.py` | `c650a61a37c8ef2b1c9817ce598babe4` | 4998 |
| `step2_run_bayesprism_tcga.R` | `00dd0618941d86acbeb3a5c42e4339a7` | 4873 |

Those two are unchanged from the deposited originals. The three below were
repointed on 2026-09-03; their pre-migration bytes are the md5 in the second
column, and are on disk beside them as `*.bak_20260903`.

| file | md5 as imported | changed |
|---|---|---|
| `step3_survival_by_fraction.py` | `5442deea407c4ea3691ec2158d3f9f6d` | input directory, output directory |
| `step4_extract_epi_cd8.R` | `f4d58394dd588742a3f3f3002d914402` | input directory, `setwd` |
| `survival_nfkb_cytokines.py` | `c34385dee75249a8459d7ea923bb7222` | input directory, `TCGA_BASE` |

Also here, and not part of the deposited set: `step1_prepare_tcga_bulk_pinned.py`
and `step4_extract_epi_cd8_pinned.R`, the two determinism fixes VERDICT.md
records. They used to exist only inside `rerun/` and `rerun_step1/`, which
`run_cpu_pipelines.sh`'s `stage()` deletes on its way in, so each was one
pipeline run away from being lost. They live here now and are staged from here.

## Migration of 2026-09-03

Author's ruling: migrate the chain off the temp workspace. Check 11 of
`03_Revised_Panels/verify_panel_provenance.py` reported eleven references
reading live inputs out of `submission-tree/temp-workspace/02162026_BayesPrism`,
a directory that can be cleared at any time and takes Figure 3B and 3E's input
path with it when it is.

Nine data files were **copied** - never moved - into the pipeline's declared
home, `paths.py` `TCGA_BAYESPRISM_DIR` =
`submission-tree/02_Preparation_for_Panels/BayesPrism_TCGA`, which already held two of
them. Every md5 was taken at the source and again at the destination and
matched, and every one agrees with the 31 Aug 2026 table in `VERDICT.md`. The
manifest is at the end of that directory's `PROVENANCE.txt`. The temp workspace
was re-checked afterwards: all 33 files, every md5 identical, nothing deleted.

`run_cpu_pipelines.sh` was repointed with them. It staged from
`$TEMP/02162026_BayesPrism`, and check 11 could not see that: its `ABS_PATH`
regex only matches a literal `/path/to/machine/...` path and this one was built out
of a shell variable. Leaving it would have made the migration cosmetic.

**Nothing computed moved.** Proof:

- `step4_extract_epi_cd8.R`, re-run against the migrated posterior, reproduces
  both deposited tables byte for byte -
  `tcga_bayesprism_epithelial_expression.tsv` `2c81930ccb80241ea8407177d19558da`
  and `tcga_bayesprism_cd8_expression.tsv` `6a16242fc9a1cc77bca8010c00f97add`.
  That is the artefact Figure 3B and 3E are drawn from.
- `step3_survival_by_fraction.py` and `survival_nfkb_cytokines.py` were run
  twice from copies differing in one line - the input directory - and produced
  identical stdout and identical PNG md5s.
- All 441 panel artefacts fingerprinted by `10_Reproduction/capture_baseline.py`
  are unchanged (`objtable-1`), Figure 3B and 3E re-run included.
- Each comparison was positive-controlled: perturbing one cell of
  `tcga_bayesprism_fractions.tsv` across the median moved the printed means,
  two log-rank P values and two of the three PNG md5s, and shifting one path
  coordinate by 1.5 units in `03_C/tcga_ceacam_cd274.svg` made
  `verify_equivalence.py` name that file and report NOT EQUIVALENT.

## Two pre-existing faults found while doing this, neither caused by it

**1. `step3_survival_by_fraction.py` crashed partway through, and did before.
Repaired 2026-09-03 on the author's ruling.**

Sections 1-3 (fibroblast, MoMac and combined fraction KM) completed. Section 4
raised `ValueError: Values must be numeric` from `lifelines`: `merged` carries
its own row index after the `.merge()` in `load_data()`, and
`merged_fib.join(fib_il6, how='inner')` joins on that index, not on the sample
ids `fib_il6` is indexed by. Nothing matched, the frame came back empty, and
sections 4 and 5 never ran. The traceback was byte-identical before and after
the migration, so the defect is in the deposited script, not in the repoint.

The repair is two lines, both joins, and nothing else: `submitter_id` is already
a column of `merged`, so the join is given it -
`join(fib_il6, on='submitter_id', how='inner')`. Index alignment only. Same
table, same gene, same log2, same median split. The pre-repair bytes are on
disk as `*.bak2_20260903`.

**There is no deposited artefact to check the repair against, for any section.**
`find -L` over every tree of the project finds **zero** copies of all five
files this script writes - `survival_fibroblast_fraction.png`,
`survival_momac_fraction.png`, `survival_fib_momac_combined.png`,
`survival_fib_il6.png` and `survival_momac_il1b.png`. None of them has ever
existed on disk. An earlier version of this file said the last two were
"produced by something that did work"; that was wrong, and is corrected here.
The `survival_*.png` files that ARE in the deposited working directory
(`survival_individual_*.png`, `survival_momac_nfkb_4gene.png`,
`survival_momac_nfkb_5gene.png`) are written by no script in `scripts/`.

So the output below is **new**. It is not a reproduction of anything, and a
future reader must not read it as one. Nothing printed in the paper depends on
it: `03_Revised_Panels/PROVENANCE.csv` has no row for any of the five files, and
the only references to them anywhere in the three round trees are this script,
its `rerun/` copy, their backups, and this file.

What the repaired script produces, 3 Sep 2026, from the migrated inputs:

| section | split | n high | n low | log-rank P |
|---|---|---|---|---|
| 1 Fibroblast fraction | median | 183 | 183 | 0.1142 |
| 2 MoMac fraction | median | 183 | 183 | 0.8671 |
| 3 Fib + MoMac combined | median | 183 | 183 | 0.0641 |
| 4 Fibroblast-specific *IL6* | median of log2(x+1) | 184 | 182 | **0.0567** |
| 5 MoMac-specific *IL1B* | median of log2(x+1) | 204 | 162 | **0.2541** |

Sections 1-3 are byte-identical to what the script produced before the repair -
all three PNG md5s and all 29 lines of stdout. Sections 4 and 5 are new.

The uneven splits in 4 and 5 are ties, not a fault: the script's rule is
`>= median`, and 2 samples sit exactly at the IL6 median and 46 at the IL1B
median, so they all go to the high group. Both joins now match 366 of 366
analysed samples. Both P values were reproduced independently, outside the
script, straight from the tables. A positive control - moving one sample's IL6
across the median, in a sample that is in the analysed set - moved section 4's
P from 0.0567 to 0.0616 and changed only `survival_fib_il6.png`, leaving the
other four unchanged.

**2. `survival_nfkb_cytokines.py` could not run at all.** `TCGA_BASE` reached
the clinical table through `upstream-pipeline/.../98_External/Bulk/02_TCGA_STAD`, a
symlink into `submission-tree/01_Raw_Inputs/04_External`, and `04_External` was renamed
`02_External`. This is fault 3 of `VERDICT.md`, already recorded there for
`step1_prepare_tcga_bulk.py`. It is repointed here at the same file under the
name it now has, exactly as `step1_prepare_tcga_bulk_pinned.py` repoints it -
same bytes, different route. With that fixed the script runs, and both sides of
the input-directory comparison above used the fixed route, so the comparison
stays single-variable.

## Correction: the panel letters in the first version of this file were wrong

It said "Figure 3C, 3D", carried over from the `PIPELINES` table in
`import_pipelines.py`. Those are build directory letters. The project's second standing rule:
a panel directory's letter is not the printed panel letter, and
`03_Revised_Panels/PROVENANCE.csv` is the authority. Build directory `03_C`
holds printed panel **B** and `03_D` holds printed panel **E**. Printed panels
3C and 3D are `03_E/create_cd8_umap.py` and `03_B/create_ceacam_tex_scatter.py`,
which read the h5ads and have nothing to do with this pipeline.

Copied verbatim except where this file says otherwise; any further change made
for determinism is recorded in `run_cpu_pipelines.sh`, never by silently
editing these.
