# Supplementary figures S7–S11 — this tree IS the build path

The opposite of `../Main_Figures/`. These five figures are new in this revision and are
genuinely rebuilt from the scripts below:

    S*/S*_*/create_S*_*.py
      -> S*/S*_*/*.svg
      -> assemble_new_supplementaries.py
      -> _assembled/S*.pdf                           shipped by 06_Submission_Package

Editing a panel script here **does** change what ships, once you re-run the script and
then `assemble_new_supplementaries.py`.

## Here the directory letters do match

`S9_C` holds panel C of Figure S9. That is true throughout S7–S11, and it is *not* true
of the main figures, where twenty-five printed panels sit in a differently named
directory. Do not carry the habit across. `../PROVENANCE.csv` covers both.

## The other supplementary figures

- **S1–S6** are carried over from `../../00_GROUND_TRUTH/figures/`. S1 is the one
  exception, and **no panel of it is regenerated**: `../Supplementary_Fixes/patch_S1_sample_labels.py`
  replaces the sample-axis tick text of panels **B and C** in the submitted PDF and
  changes nothing else in the file. Those two panels are where the submitted figure
  printed internal specimen identifiers. It raises rather than falling back if any
  specimen is unmapped — do not add a fallback.

  `../Supplementary_Fixes/S1_H/` is **not** that fix, and the reason once given for it
  here has been retracted. It was built on the premise that panel H carried the raw
  identifiers; `00_Data_Audit/FINDINGS.md` §12.13 (31 Aug 2026) established that they
  were in panel C and that panel H's labels were already study IDs. The script is
  correct and reproduces byte-identically, but nothing places its `panel_S1_H.svg` into
  any figure, and panel H ships carried over like the rest of S1. Its only consumer
  anywhere is `assemble_S1.py` inside the two **generated** release trees, which does
  not reproduce the submitted page. See `../SUPPLEMENTARY_AUDIT.md` fault 4, and the
  S1 rows of `../PROVENANCE.csv`, which already record this.
- The analyses behind S7–S11 live in `../../02_New_Analyses/`, one numbered module per
  reviewer point. Several panel scripts here are thin wrappers that read a module's
  `outputs/*.csv`; the module is the place to change a number, not the panel.

## Reading expression data

From `.raw`, never `.X`. Eight input h5ads carry a double normalisation from 2026-07-30
that left whole cell rows NaN — `00_Data_Audit/FINDINGS.md`, sections 1 and 7. Scripts
here raise if `.raw` is absent rather than falling back.
