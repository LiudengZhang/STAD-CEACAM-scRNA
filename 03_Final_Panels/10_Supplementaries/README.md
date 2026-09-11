# Supplementary figures S7–S9 — this tree IS the build path

The opposite of `../Main_Figures/`. These three figures are new in this revision and are
genuinely rebuilt from the drivers below:

    _drivers/draw_*.py
      -> S*/S*_*/*.svg
      -> assemble_new_supplementaries.py
      -> _assembled/S*.pdf                           shipped by 06_Submission_Package

Editing a driver here **does** change what ships, once you re-run it and then
`assemble_new_supplementaries.py`.

The supplement is S1–S9. The eleven panels of this tree are:

    S7  CEACAM5 versus CEACAM6            A, B, C
    S8  Spatial confounders               A, B
    S9  MoMac identity and NF-κB          A, B, C, D, E

One further panel is drawn here and does not belong to a figure of this tree:
`_drivers/draw_adaptive_immune_resource.py` draws **S3 E**, which
`../Supplementary_Fixes/patch_S3_add_adaptive.py` composes below the submitted S3 page.

## A panel is drawn at the size it prints at

Every driver draws its panel at its millimetre print size and the assembler places it
1:1; `check_scale()` refuses to assemble if any placement is not 1:1. Type set at 7 pt
prints at 7 pt. Do not reintroduce a scale factor at draw time and do not let the
assembler fit a panel into a box.

A driver imports the analysis module that owns its numbers, reads the tables that
module's `main()` already wrote, and calls its plotting entry points. It recomputes no
statistic. Run the analyses before the drivers.

## Here the directory letters do match

`S9_C` holds panel C of Figure S9. That is true throughout S7–S9, and it is *not* true
of the main figures, where twenty-five printed panels sit in a differently named
directory. Do not carry the habit across. `../PROVENANCE.csv` covers both.

## The other supplementary figures

- **S1–S6** are carried over from `../../00_GROUND_TRUTH/figures/`, three of them with a
  patch applied in place by `../Supplementary_Fixes/`:
  - **S1** — `patch_S1_sample_labels.py` replaces the sample-axis tick text of panels
    **B and C** in the submitted PDF and changes nothing else in the file. Those two
    panels are where the submitted figure printed internal specimen identifiers. It
    raises rather than falling back if any specimen is unmapped — do not add a fallback.
    **No panel of S1 is regenerated.**
  - **S2** — re-exported with two-sided P values.
  - **S3** — `patch_S3_add_adaptive.py` draws the submitted page onto a taller canvas
    exactly as it stands and composes panel **E** below it. Nothing above the join is
    re-rendered, and `gate()` proves that pixel row by pixel row.

  `../Supplementary_Fixes/S1_H/` is **not** the S1 fix. It was built on the premise that
  panel H carried the raw identifiers; `00_Data_Audit/FINDINGS.md` §12.13 established that
  they were in panel C and that panel H's labels were already study IDs. The script is correct and reproduces
  byte-identically, but nothing places its `panel_S1_H.svg` into any figure, and panel H
  ships carried over like the rest of S1. Its only consumer anywhere is `assemble_S1.py`
  inside the two **generated** release trees, which does not reproduce the submitted
  page. See `../SUPPLEMENTARY_AUDIT.md` fault 4, and the S1 rows of `../PROVENANCE.csv`,
  which already record this.
- The analyses behind S7–S9 live in `../../02_New_Analyses/`, one numbered module per
  reviewer point. The module is the place to change a number, not the driver.

## `_restyled/`

Holds the cnsplots restyle of the previous twenty-six-panel S7–S11 set. The tree it
restyles no longer exists and the restyle is ruled out of the code release, so nothing
it draws is printed. It is kept only because the release tooling still addresses it by
path; it is archived under `07_Archive/`, and it should be removed in the same step that
regenerates `RELEASE_UNIVERSE.csv` and `RELEASE_GAPS.csv`.

## Reading expression data

From `.raw`, never `.X`. Eight input h5ads carry a double normalisation that left whole
cell rows NaN — `00_Data_Audit/FINDINGS.md`, sections 1 and 7. Scripts here raise if
`.raw` is absent rather than falling back.
