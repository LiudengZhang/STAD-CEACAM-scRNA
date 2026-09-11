# Figure 5 panel F — superseded, not broken

Ruling-date: 2026-09-09
Ruled-by: the author (CIR-26-0753-ET major revision)
Panel: Figure 5 **F** — the printed letter, looked up in
`03_Revised_Panels/PROVENANCE.csv`. The directory is `05_H`; do not read the
directory as the letter.
Script: `create_cytokine_dotplot.py`
Marker: this file, **not** `KNOWN_BROKEN.md`. See `RULES.md` rule 1.

## What was superseded

The dot **colours** of the panel as printed in the submitted figure
(`00_GROUND_TRUTH/figures/Figure 5.pdf`). The redrawn panel does not reproduce
them and is not meant to. Nothing else moved: same four genes, same thirteen
cell types, same row order, same `major_cell_type` filter, same `cmap='Reds'`,
same `standard_scale='var'`, same strings.

## Why

The printed colours were computed from a damaged matrix.

`full_dataset.h5ad` has no `.raw` — `06_Clean_Data/build_clean_h5ad.py` promoted
`.raw.X` to `.X` — and for this file that `.raw` was the matrix normalised a
second time on 2026-07-30, which `00_Data_Audit/FINDINGS.md` 12.11–12.12
measures and reproduces. The panel drew through `use_raw=True`, so its colours
came from that matrix.

That is a measurement, not an inference: the 52 dot fills of the previous SVG
were read back out of the file, the `Reds` ramp inverted, and they reproduce the
group means the damaged matrix gives.

## What it was redrawn from

`layers['counts']` of the deposited object — intact integer UMI counts,
`normalize_total(target_sum=1e4)` then `log1p`. That transform is the project's
own, read out of FINDINGS 12.12 rather than assumed from scanpy's default, and
validated before use: redoing it from the counts layer of
`Neutrophils_sound.h5ad` reproduces that file's own `.X` to max |diff| 0.000000.

The input is `FULL_DATASET_DEPOSIT_H5AD` (`00_Config/paths.py`), not
`FULL_DATASET_H5AD`, which resolves into the submission-tree input tree to a file with
no counts layer. The two matrices were measured identical over all 542,121 cells
and all four plotted genes — max |diff| 0.0000000000, 0 cells differing,
`var_names` identical.

## What moved, exactly

    dot size    percentage of expressing cells. Depends only on the zero
                pattern, and a monotone transform cannot move it. Measured:
                max |delta pct| = 0.0000000000 over 13 groups x 4 genes.

    colour      each group's mean, scaled per gene by standard_scale='var'.
                Moves by at most 0.069 (TNF), 0.121 (IL6), 0.090 (IL1B),
                0.094 (IL1A). IL6's highest group becomes Fibroblast where the
                shipped panel printed B cells. The highest group for TNF (DC)
                and for IL1B and IL1A (both MoMac) is unchanged.

## Why there is no Retest-by date

There is nothing to retest. The panel will not reproduce the printed one again,
by design and under the ruling above. A `Retest-by:` line here would file an
authorised correction as a defect. The two senses of `reproduces_published = no`
are set out in `RULES.md` rule 1 and enforced by
`03_Revised_Panels/verify_panel_provenance.py` checks 3 and 10.

## Record

- Previous version of the panel and its script:
  `07_Archive/2026-09-09_fig5F_sound_redraw/`
- The deposited object as it stood before its `.X` was recomputed from counts:
  `07_Archive/2026-09-10_before_X_recompute_from_counts/`
- `PROVENANCE.csv`, Figure 5 panel F: `reproduces_published = no`, with this
  mechanism in its note.
