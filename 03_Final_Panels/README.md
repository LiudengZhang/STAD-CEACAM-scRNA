# Main figures 1–6 — the panel sources, and how the shipped pages are built

Until 2026-09-09 this tree was archive and the figures in the paper were the submitted
PDFs edited in place (`../patch_figure_annotations.py` → `_patched/`). That path is
retired. Since 2026-09-11 (Figures 2–5) and 2026-09-15 (Figure 1) **every panel of
Figures 1–5 is redrawn here at the size it prints at**, and the pages are assembled
from the redrawn panels:

    0X_Figure_X/<dir>/create_*.py            one script per panel, 1:1, on
                                             00_Config/panel_style_cns (cnsplots),
                                             numbers from data/*.csv (cnsfig.cache)
      -> ../../12_Figure_Refactor/rebuild_panels.py N --only <dir>
      -> ../assemble_slotted.py N            places each drawing into its slot from
                                             ../panel_rects_v2.csv / slot_subrects_v2.csv
                                             (written by 12_Figure_Refactor/build_grid_v2.py)
      -> _slotted/Figure_N.{svg,pdf,png}
      -> ../build_shipped_figures.py         collects the pages
      -> _shipped/Figure_N.pdf               shipped by 06_Submission_Package

Every page is 171.10 × 229.31 mm. Figure 1 was submitted as a 254 mm landscape page and
is re-paged onto the same width (`build_grid_v2.REPAGED`): A, the study-overview
artwork from `_panel_1A/build_figure.py`, across the top; B (UMAP) and C (composition,
upright bars, key below) on one row at one height, 62 mm, since 2026-09-16. Figure 6 is
supplied as submitted (`_supplied/Figure_6.pdf`).

**A redraw changes the drawing only.** Not one number, gene set, cell selection,
statistic or string moves without a declaration in `00_Config/shared/labels.py`;
`../../10_Reproduction/compare_panel_content.py` is the proof, run for every redrawn
panel. The departures the author ruled on (one box style with no individual points,
2026-09-16; abbreviated labels, 2026-09-10; and others) are recorded in
`../PROVENANCE.csv` row by row and in `../../RULES.md`.

`assemble_figure_*.py` beside the panel directories are the **pre-submission**
assemblers. They are not on the build path and must not be used to rebuild the paper:
`assemble_figure_2.py` emits panels A–Q where the paper prints A–N, and
`assemble_figure_1.py` prints the submitted landscape page. They are kept as the record
of how the submitted pages were made and travel with the code release for that reason.

## The directory letter is not the printed panel letter

Twenty-five printed panels live in a directory whose name says something else:
Figure 5 panel H is in `05_F`, panel J is in `05_I`, panel G is in `05_C`; Figure 3
panel B is in `03_C`; Figure 2 panel G is in `02_H`. Three different Figure 3 panels
share the filename `ceacam_spatial_boxplot.png`.

**Look every panel up in `../PROVENANCE.csv`.** Do not infer one from a directory name.

## Figure 5A: the printed panel is correct

`05_Figure_5/05_A/` reproduces the published panel against the GSEA table the panel was
made from, `MoMac_mast_prerank_gsea.csv`: the top nine hallmarks by |NES| are the printed
nine, in the printed order, matching the published bars to 0.00014 NES. It reproduces
nothing at all if the script is pointed at a different run — `GSEA/post/`, built on the
doubly normalised matrix — which is what an input path moved into an archive will do to
it. When a panel will not reproduce, suspect its input path before the figure;
`../verify_panel_provenance.py` check 9 is what looks at where every panel script reads
from.

## Rebuilding a panel

Edit its `create_*.py`, run `rebuild_panels.py N --only <dir>` (it refuses a script that
exits 0 without writing a drawing), `assemble_slotted.py N`, `build_shipped_figures.py`,
then the gates in `../../RULES.md` "Before shipping" — `sweep_panels.py`,
`sweep_pages.py`, `check_figure_annotations.py`, `verify_panel_provenance.py` — and the
comparator against the previous drawing. Update the panel's row in `../PROVENANCE.csv`
with a dated sentence saying what changed.
