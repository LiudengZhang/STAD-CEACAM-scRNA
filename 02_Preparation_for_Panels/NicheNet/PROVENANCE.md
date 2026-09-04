# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# NicheNet — where this code came from

> # ⛔ DO NOT "FIX" THE `colnames` LINE. READ THIS FIRST. ⛔
>
> **File:** `Round_5/02_Preparation_for_Panels/NicheNet_v2_postNR/02_Scripts/02_run_nichenet_c3mac.R`
> **Line 116:** `nichenet_targets <- colnames(ligand_target_matrix)`
>
> **This line is wrong on its face and is being kept anyway. The author ruled on
> 2026-09-03: keep the code exactly as it is. Figure 5 panel G does not change.
> DO NOT FIX IT.**
>
> If you have just noticed this defect and are about to correct it: stop. It has
> already been found, reviewed, and deliberately retained. Correcting it would
> silently move a published panel. This project has been bitten by exactly that
> failure mode twice, both times on Figure 5, both times retracted.

## The defect, precisely

The prior model `ligand_target_matrix.rds` is oriented **25,345 target genes
(rows) × 688 ligands (columns)** — the standard nichenetr orientation. Line 116
therefore assigns the **688 ligand symbols** to a variable named
`nichenet_targets`, not the 25,345 target genes.

Measured in R on 2026-09-03, not inferred:

    dim(m)                          25345 x 688
    head(rownames(m), 8)            A1BG, A1BG-AS1, A1CF, A2M, A2M-AS1, A2ML1, A2MP1, A3GALT2
    head(colnames(m), 8)            CXCL1, CXCL2, CXCL3, CXCL5, PPBP, CXCL6, CXCL8, CXCL9
    intersect(colnames(m), lr$from) 688 of 688   <- every column is a ligand
    intersect(rownames(m), lr$from) 688 of 25345

The rows are the alphabetical genome-wide gene space; the columns are the
ligands, and all 688 of them appear in `lr_network$from`. The correct expression
for the target space would be `rownames(ligand_target_matrix)`.

Consequences, both retained:

1. **The gene set of interest is drawn from the ligand column space, not the
   target space.** Lines 116–122: `receiver_nichenet_genes` is
   `intersect(colnames(receiver_expr), nichenet_targets)`, and the top-250
   most variable of those become `target_genes`, which is what
   `predict_ligand_activities(geneset = ...)` receives at line 137. The geneset
   is therefore up to 250 of the 688 ligands, ranked by variance in the
   receiver, rather than 250 receiver target genes.
2. **`ligand_targets.csv` is a ligand × ligand submatrix.** Line 177,
   `ligand_target_matrix[top_ligands, target_genes]`, indexes rows by ligand
   symbol against a matrix whose rows are target genes, and columns by
   `target_genes` (ligand symbols) against columns that are ligands. The file
   written at line 186 is consequently not a ligand→target table. Visible in the
   deposited files: `results/B_cells/ligand_targets.csv` has `IL10` and `IL15`
   in its `target` column, which are ligands.

Line 78 (`intersect(expressed_genes_sender, rownames(ligand_target_matrix))`)
has the mirror-image slip — it intersects against the 25,345 target genes where
the ligands are wanted — but it is **harmless**, because line 82 immediately
re-intersects with `unique(lr_network$from)`, which restricts the result to the
688 ligands regardless. That is why `expressed_ligands` comes out as a sensible
106. It is recorded here only so a future reader does not "discover" it
separately and conclude the ligand set is wrong too. It is not.

`predict_ligand_activities` itself receives the matrix in its correct standard
orientation and uses it correctly; the defect is confined to lines 116 and 177.

Note also that the console message at lines 50–51 prints the two dimensions
with the labels transposed ("ligands x targets"), which is how the orientation
was mis-read in the first place.

## The ruling

**Date: 2026-09-03. Author's decision: keep the code exactly as it is. Do not
fix line 116, and do not fix line 177.**

`c3mac_ligand_activities_raw.csv` — the deposited table that
`00_Config/paths.py:LIGAND_ACTIVITIES_RAW` hands to the panel — was produced by
this behaviour. **Printed Figure 5 panel G rests on it** (`03_Revised_Panels/
PROVENANCE.csv` row 51: figure 5, printed panel G, build directory `05_C`,
script `05_C/create_panel_c_ligand_violin.py`, `reproduces_published = yes`,
adjudicated 2026-09-01). Changing line 116 changes the gene set, which changes
every AUROC/AUPR/Pearson in that table, which changes the top-15 ligands the
violin panel draws. That is a published number moving. It is not permitted here.

The defect was **reviewed and deliberately retained**, not overlooked. This is a
record of a decision, not a bug report awaiting action. Anyone who wants it
changed needs the author's word, in writing, first — and a redraw of Figure 5G
is explicitly *not* that: a restyle is a visualisation-only change and may not
move a value (see the project `CLAUDE.md`, rule 3).

## The manuscript is being corrected instead

The fix is in the prose, not the code. The Methods text is being amended
separately, concurrently, in `04_Manuscript_R1/01_Main_Text/edits.py`, so that
it describes the gene set the code actually builds. Do not edit that file from
here.

## Where the real pipeline lives

The canonical run is **not** the Round_4 prototype this directory's `scripts/`
copy was taken from — that prototype never produced any `c3mac_*` output. It is:

    Round_5/02_Preparation_for_Panels/NicheNet_v2_postNR/
      01_Config/config.yaml              both scripts read this; no CLI arguments
      02_Scripts/01_prepare_data.py      conda run -n Liudeng_Python_310
      02_Scripts/02_run_nichenet_c3mac.R conda run -n r_demo
      PROVENANCE.txt                     dated 2026-02-22
      results/c3mac_ligand_activities_raw.csv
                                         md5 d6824aa593a2c64e964373d29208e6ee,
                                         identical to the deposited copy in
                                         Round_5/02_Preparation_for_Panels/NicheNet/

A full parameter recovery with file:line evidence is at
`/path/to/home/graphst_repro/nichenet_recovery/recovered_parameters.yaml`.

**"NicheNet_v2" in these directory names means the project's *second* NicheNet
analysis. It does not mean the NicheNet v2 model.** The prior actually used is
the **v1 human** model, `ligand_target_matrix.rds`, 25,345 × 688,
md5 `f79a1b2a7bf4c401319450e813ce62fe`. The v2 `..._nsga2r_final.rds` prior was
tried and rejected — it gave LILRB4/TFF1 as top ligands
(`Round_5/.claude/plans/completed/rerun_nichenet_c3mac_sender.md`).

## Where the prior models live - migrated 2026-09-03

They used to be read out of an archive.
`Round_5/02_Preparation_for_Panels/_archived_NicheNet_v2/00_Databases/` was
named by lines 37-39 of **both** live copies of `config.yaml` - the canonical
one at `NicheNet_v2_postNR/01_Config/` and this directory's
`scripts/01_Config/` copy, which is the one `02_Upstream/run_cpu_pipelines.sh`
actually stages. Check 11 of `03_Revised_Panels/verify_panel_provenance.py`
failed six times on it, three lines by two copies. That is the Figure 5A
failure shape exactly: a live input reached through an archive path.

They now live at

    Round_5/02_Preparation_for_Panels/NicheNet/00_Databases/
      ligand_target_matrix.rds   f79a1b2a7bf4c401319450e813ce62fe  123,235,892 B
      lr_network.rds             0dcac9202cdbfb091a6628d31ff2efed       26,250 B
      weighted_networks.rds      eedca2cb2a286c557f0b992296c2d8ef   35,326,653 B

declared as `00_Config/paths.py:NICHENET_DB_DIR` - which is what the config's
own comment at line 36 had promised since February and which did not exist
until now. They were **copied**, never moved: the archive is byte-for-byte as
it was, all 91 files re-md5'd afterwards, nothing deleted. Only data was put
there; the scripts stay where they are, because `update_release.py` syncs
`*.py`/`*.R`/`*.yaml` out of `02_Preparation_for_Panels` into the published
repository and these are pipeline code, not preparation data.

**Nothing computed moved.** `02_run_nichenet_c3mac.R` - byte-identical to the
deposit, md5 `8ea382559eec432bd53a8f89650bb951`, line 116 untouched - was
re-run from the repointed config against the deposited `prepared_data/`, so
the three database paths were the only variable. All 42 outputs are
md5-identical to `NicheNet_v2_postNR/results/`, and Figure 5G redraws to a
byte-identical PNG. `VERDICT.md` and
`Round_5/02_Preparation_for_Panels/NicheNet/PROVENANCE.txt` carry the digests
and the positive controls.

`weighted_networks.rds` is declared by `config.yaml` but read by neither
script - `02_run_nichenet_c3mac.R` loads only `ligand_target_matrix` and
`lr_network`. It was migrated with the other two so that the declared input set
is complete and nothing later reaches back into the archive for it.

One reference was deliberately **not** changed:
`scripts/01_Prototype/run_nichenet_analysis.R:54` hard-codes a `db_dir` under
`Round_4/02_Playground/02_NicheNet/00_Databases`, a path that no longer exists
(the live one is `.../02_Playground/_archived/02_NicheNet/00_Databases`, and it
holds the *nsga2r* prior, not this one). That file is the dead Round_4
prototype kept verbatim as history; it produced none of the deposited output
and check 11 does not flag it. Editing it would falsify the record of what the
prototype was.

Produces: ligand activities for the C3_Mac_Inflam_IL1B sender across 13 receivers.

Consumed by: **printed Figure 5 panel G** (build directory `05_C`). Earlier
copies of this file said "Figure 6"; that was wrong. Never infer a panel letter
from a directory name — look it up in `03_Revised_Panels/PROVENANCE.csv`.

## The archived Round_4 prototype

`scripts/01_Prototype/` below is a verbatim copy of the dead Round_4 prototype
`Round_4/01_Round_4.2_Standardized_Pipeline/03_Cell_Cell_Interaction/02_NicheNet_Analysis`.
It is kept as history. It is **not** what produced the deposited artefacts.

| file | md5 | bytes |
|---|---|---|
| `01_Prototype/core_nichenet_pipeline.py` | `a28011b26d8e6540a252818cb50f16b5` | 5017 |
| `01_Prototype/data_preparation.py` | `f2fc761a3f469aeca5c9b1688c3de164` | 8430 |
| `01_Prototype/nichenet_runner.py` | `c31abfbfb375fd852abe81be1a5cbf43` | 3463 |
| `01_Prototype/utils.py` | `f4224dc5f2323cdcbaab089247e00d58` | 5186 |
| `01_Prototype/visualization.py` | `bf8f3c5587787f830d60fd121f427d5f` | 3747 |
| `01_Prototype/run_nichenet_analysis.R` | `56e0cc7060390d12312dff483693a6e2` | 6200 |

Copied verbatim; any change made for determinism is recorded in `run.sh`, never
by editing these.
