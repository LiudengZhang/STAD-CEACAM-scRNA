> **Upstream record.** Paths in this document refer to the upstream
> processing pipeline that produced the deposited intermediates; it is not
> part of this release's run path. The files under `upstream/` are included
> as a record of how the inputs were produced and are not executed by
> `./run` or `_run_all_panels.sh`.

# PRJEB25780 (TIGER): the two tables that are NOT in the Zenodo archives

Added 2026-10-05.

Three shipped panels, the two-sided sweep and the metaprogram external
validation read two PRJEB25780 tables. Neither is in any Zenodo archive and
neither may be, because neither is ours to redistribute. Both are regenerated
here, exactly, from TIGER's own two downloads.

## What a reader needs, and where it comes from

TIGER — <http://tiger.canceromics.org/> — publishes the PRJEB25780 cohort as
two R objects. Download both for dataset `STAD-PRJEB25780` and put them in a
directory of your choosing as `00_Input/`:

| file | bytes | md5 |
|---|---|---|
| `STAD-PRJEB25780.Response.Rds` (expression) | 16,206,715 | `330680e4f3ea62b208de6c90156f3284` |
| `STAD-PRJEB25780.Response-2.Rds` (sample metadata) | 1,036 | `f6aecc32811b6cbdf27a3a380016f537` |

The digests are the point of the table. TIGER is a curated database, not an
archive with versioned releases, so "the PRJEB25780 metadata" is not by itself
a statement about which bytes you hold. These two are the bytes the paper used.

## Then run

```bash
export TIGER_BASE_DIR=/path/to/that/directory      # it holds 00_Input/
Rscript scripts/01_prepare_expression_data.R
Rscript scripts/02_run_estimate.R
```

Both write into `$TIGER_BASE_DIR/02_Tumor_Purity_Analysis/01_ESTIMATE_Results/`
and log into `$TIGER_BASE_DIR/03_Logs/`. Needs R with `data.table`, `dplyr`,
`tidyr` and `estimate`; 33 seconds for both on one core.

## What comes out, and where the code expects it

| produced | md5 | the code reads it as |
|---|---|---|
| `expression_matrix_for_estimate.txt` | `f1782308b7f240f178e841dcea0a174d` | `02_External/Bulk/PRJEB25780_Tiger/expression_matrix.txt` |
| `metadata_with_purity_scores.tsv` | `2e012a4d6294fb424ff1caaffd1ac044` | `00_Config/paths.py:TIGER_META`, i.e. `02_External/Bulk/PRJEB25780_Tiger/metadata_with_purity.tsv` |

Copy or link them to those two names under `$STAD_RAW_INPUTS` and the panels
and analyses that read them run unchanged.

**Reproduced byte for byte on 2026-10-05**, with the copies in `scripts/`
here, into a directory that had never held an output: both md5s above came
back identical. Re-measured after the byline was removed, because a proof about
a file that has since changed is not a proof about the file that ships. That measurement is why the two tables are not in the deposit.
The question "the ESTIMATE columns are ours, so is the joined table ours to
publish?" never has to be answered — the file regenerates, so nothing is lost
by leaving it out.

## What the ESTIMATE columns are, and the record that said otherwise

`metadata_with_purity.tsv` is TIGER's sample metadata with four columns added:
`StromalScore`, `ImmuneScore`, `ESTIMATEScore` and `TumorPurity`. They are
**computed here**, by `02_run_estimate.R`, from TIGER's expression matrix
through the `estimate` package. Two records in the release said otherwise and
were corrected the same day:

- the capsule README called the table TIGER metadata "joined to that record's
  ESTIMATE scores", which reads as though TIGER published the scores;
- `External_Public_Inputs/PROVENANCE.txt` said "no script in the project
  produces them".

Both were settled by digest, not by argument. The deposited
`expression_matrix.txt` is byte-identical to this pipeline's
`expression_matrix_for_estimate.txt`, and the deposited
`metadata_with_purity.tsv` to its `metadata_with_purity_scores.tsv`.

## The two lines that were changed

Each script carried a single hardcoded `BASE_DIR` naming
`upstream-pipeline/98_standardized-pipeline/...` — a directory that no longer
exists, since that tree was renamed `01_standardized-pipeline`. So
both scripts, as kept, could not run anywhere. `BASE_DIR` now reads
`TIGER_BASE_DIR` from the environment and stops if it is unset.

The second is the original's `# Author:` header line, which credited an AI
assistant. `scan_publish_safety.py` refuses an assistant citation in a
deposited file, and it caught both scripts on the first scan after they were
added - and then caught them again, because the comment written to explain the
removal quoted the line it had removed. A citation is a citation either way.
The gate was right twice. Nothing else differs; the originals digest

    01_prepare_expression_data.R   987626f7abb80f2e21ad2c646e62740a
    02_run_estimate.R              46fbf78542dca5d46c35318f90ff2bf5

and are at `upstream-pipeline/01_standardized-pipeline/00_Set_Up/`
`00_Database_and_External_Datasets/98_External/Bulk/01_PRJEB25780_Tiger/`
`01_Scripts/` in the working tree.

## Not deposited, and not a gap

`08_Deposit_Staging/build_external_inputs.py` refuses to put either table into
`External_Public_Inputs.zip`, by content as well as by name, and carries a
control that plants the forbidden content under an innocent name. A licence
condition that depends on nobody re-adding a file by hand is not a condition.
