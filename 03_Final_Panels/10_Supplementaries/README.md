# Supplementary figures S1–S10 — this tree IS the build path

Every supplementary figure the revision ships is drawn and assembled here:

    S<n>_*/S<n>_<letter>/create_*.py     S1–S7: one standalone script per panel
    S<n>_*/create_*.py                   (S6 and S7: one script at the figure level)
    _drivers/draw_*.py                   S8–S10 and S4 E: drivers over the analyses in
                                         ../../02_New_Analyses/
      -> S<n>_*/S<n>_<letter>/*.svg      one drawing per panel directory, at print size
      -> assemble_new_supplementaries.py rows of panels, letters from ../PROVENANCE.csv
      -> _assembled/S<n>_*.pdf           one page each, 171.10 mm wide, no taller than a
                                         main-figure page (229.31 mm); shipped by
                                         06_Submission_Package

Editing a script here **does** change what ships, once you re-run it and then
`assemble_new_supplementaries.py`.

## The numbering, since 2026-09-16

The submitted supplement was S1–S9. On the author's fifth reading the submitted S1
(QC and annotation, eight panels, 431 mm at 6 pt) was **split**: S1 is quality control
(A–E) and the new **S2 is cell-type annotation** (A–C, the former S1 F–H), and the
submitted **S2–S9 are now S3–S10**. The manuscript, the legends, the tables and the
response letter use the new numbers; the letter says so at the head of the responses.

    S1   Quality control                        A–E       (submitted S1 A–E)
    S2   Cell-type annotation                   A–C       (submitted S1 F–H)
    S3   CEACAM5/6 state and metaprogram validation  A–H  (submitted S2)
    S4   CD8+ T-cell states, adaptive composition     A–E  (submitted S3 + E)
    S5   Spatial validation                     A–F       (submitted S4)
    S6   Immune-module proportions              A–E       (submitted S5)
    S7   PD-L1 in the remaining cell types       A–I       (submitted S6)
    S8   CEACAM5 versus CEACAM6                 A–C       (submitted S7)
    S9   Spatial confounders, bulk validation   A–C       (submitted S8; B split into B, C)
    S10  MoMac identity and NF-κB               A–E       (submitted S9)

The **directories carry the printed numbers** and were renamed with `git mv` on the
same day; the caches under `data/` moved with their scripts. That is deliberate: the
assembler, `build_package.py`, `sweep_pages.py` and the release driver all read the
figure number off the directory name, and the main figures already show what a
directory that says one letter and prints another costs (`../../RULES.md` rule 2).
Here the directory letter **does** match the printed letter — `S10_C` holds panel C of
Figure S10 — and `../PROVENANCE.csv` is still the authority the assembler reads.

## A panel is drawn at the size it prints at

Every script draws its panel at its millimetre print size on `00_Config/panel_style_cns`
and the assembler places it 1:1; `check_scale()` refuses to assemble if any placement
is not 1:1. Type set at 6 pt prints at 6 pt. Do not reintroduce a scale factor at draw
time and do not let the assembler fit a panel into a box.

The S1–S7 scripts are copies of their submission-tree predecessors restyled at print size, with
their numbers in `data/*.csv` through `cnsfig.cache` (the h5ad is opened only with
`--recompute`); each was proved value-identical to its predecessor with
`../../10_Reproduction/compare_panel_content.py --figure S`. A driver imports the
analysis module that owns its numbers, reads the tables that module's `main()` already
wrote, and calls its plotting entry points; it recomputes no statistic. Run the analyses
before the drivers. The module is the place to change a number, not the driver.

Box plots throughout take one style, `00_Config/cnsfig/boxes.py` (`draw_boxes`,
`bracket`), and show no individual points — the author's ruling of 2026-09-16.

## History

- **S1–S6 until 2026-09-15 (evening)** were the submitted pages carried over from
  `../../00_GROUND_TRUTH/figures/` and patched in place by `../Supplementary_Fixes/`
  (S1's sample-axis labels, S2 D's two-sided P, S3's panel E). That tree is **retired**:
  the study IDs, the two-sided P and panel E are drawn, not patched.
- **S1 on two pages, the night of 2026-09-15 only.** The assembler's `pages` mechanism
  (`CONTINUED_NOTE`) is kept but no figure uses it.
- `_restyled/` (the cnsplots restyle of the pre-2026-09-08 S7–S11 set) is gone from this
  tree; it is archived under `07_Archive/`.

## Reading expression data

From `.raw`, never `.X`. Eight input h5ads carry a double normalisation that left whole
cell rows NaN — `00_Data_Audit/FINDINGS.md`, sections 1 and 7. Scripts here raise if
`.raw` is absent rather than falling back.
