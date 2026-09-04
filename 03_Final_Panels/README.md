# Main figures 1–6 — this tree is ARCHIVE, not the build path

**Running anything under `0X_Figure_X/` does not change what ships.** It never has.

The figures in the paper are produced by editing the submitted PDFs in place:

    ../../00_GROUND_TRUTH/figures/Figure N.pdf
      -> ../patch_figure_annotations.py
      -> _patched/Figure_N.pdf                       shipped by 06_Submission_Package

`patch_figure_annotations.py` does three kinds of thing, and nothing else:

- redacts an old annotation and redraws the exact two-sided *P* value at the same
  anchor (Figures 3 and 5, live text; Figure 2, vector outlines with no text layer);
- corrects two labels in Figure 4 (`CD16` → `FCGR3A`, `Mac` → `MoMac`);
- splices a replacement panel into a measured slot — currently Figure 1A, Figure 2D
  and Figure 5H, listed in `PANEL_SWAPS`.

The panel directories below hold the sources the pre-submission assembler consumed.
They are kept for provenance and for the code release. The assemblers
(`assemble_figure_*.py`) are **not** a way to rebuild the paper: `assemble_figure_2.py`
was re-lettered after submission and emits panels A–Q where the paper prints A–N, and
the pre-submission assembler is not on disk. Re-assembling would silently renumber the
figure and invalidate every panel reference in the text.

## The directory letter is not the printed panel letter

Twenty-five printed panels live in a directory whose name says something else:
Figure 5 panel H is in `05_F`, panel J is in `05_I`, panel G is in `05_C`; Figure 3
panel B is in `03_C`; Figure 2 panel G is in `02_H`. Three different Figure 3 panels
share the filename `ceacam_spatial_boxplot.png`.

**Look every panel up in `../PROVENANCE.csv`.** Do not infer one from a directory name.

## One panel is known not to reproduce

`05_Figure_5/05_A/` — see `05_Figure_5/05_A/KNOWN_BROKEN.md`. The printed panel is
correct; the script in that directory is not, and has twice been mistaken for evidence
that the paper was wrong.

## Rebuilding a panel

Only when a claim in the paper is wrong, never because a script disagrees with the
figure. Write the new panel to its directory, add it to `PANEL_SWAPS` with a slot
measured from the white gutters of the submitted PDF, re-run
`patch_figure_annotations.py`, and update its row in `../PROVENANCE.csv`.
