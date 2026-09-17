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

With no argument, the driver runs the main-panel workflow and the retained
pre-submission S1–S6 scripts. The `assemble_figure_2.py` through
`assemble_figure_5.py` files preserve those submitted layouts. Figure 1's
pre-submission assembler is deliberately not run because it prints the older
landscape page.

The revised main pages use the slot geometry in `panel_rects_v2.csv` and
`slot_subrects_v2.csv`. Their internal slot assembler is not included in this
release, so this repository reproduces the individual main panels rather than
assembling the revised main pages. Supplementary Figures S1–S10 are assembled
in full by `Supplementary_New/assemble_new_supplementaries.py`.

Every released panel script reads versioned data or prepared intermediates.
`PROVENANCE.csv` records its source script, input data, printed panel letter,
and reproduction verdict.
