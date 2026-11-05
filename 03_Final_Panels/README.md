# Figure panel code

This directory contains the scripts that draw the individual panels of main
Figures 1–5 and the complete builders for Supplementary Figures S1–S10.
Panel-directory letters are historical and do not always equal the printed
letter; use `PROVENANCE.csv` for the authoritative mapping.

Run one main figure's panel scripts with:

```bash
bash _run_all_panels.sh 3
```

Run the revision analyses and assemble Supplementary Figures S1–S10 with:

```bash
bash _run_all_panels.sh revision
```

With no argument, the driver runs the main-panel workflow and assembles
Supplementary Figures S1–S10 as the paper prints them. Until 2026-10-01 it ran
the superseded pre-submission S1–S6 instead; those are now `presubmission`, and
`supp` is the shipped supplement. The `assemble_figure_2.py` through
`assemble_figure_5.py` files preserve the submitted layouts. Figure 1's
pre-submission assembler is deliberately not run because it prints the older
landscape page.

This repository supplies the PANELS; the journal sets the page (the author's
ruling of 2026-10-02). The printed geometry the panels were drawn into is here as
data — `panel_rects_v2.csv`, `slot_subrects_v2.csv` and `letter_spec_v2.csv`,
one rect per printed panel, measured off the published pages
— and the five scripts that place panels onto a page and annotate it are not:
`assemble_slotted.py`, `Supplementary_New/assemble_new_supplementaries.py`,
`patch_figure_annotations.py`, `build_shipped_figures.py` and
`00_Config/shared/svg_assembler.py`. A run of `_run_all_panels.sh` produces every
panel of Figures 1–5 and of Supplementary Figures S1–S10, each at the path
`PROVENANCE.csv` declares for it.

Every released panel script reads versioned data or prepared intermediates.
`PROVENANCE.csv` records its source script, input data, printed panel letter,
and reproduction verdict.
