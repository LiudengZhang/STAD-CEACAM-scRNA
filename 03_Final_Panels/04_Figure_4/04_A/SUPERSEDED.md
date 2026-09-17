# Figure 4 panel A — superseded, not broken

Ruling-date: 2026-09-14
Ruled-by: the author (CIR-26-0753-ET major revision)
Panel: Figure 4 **A** — the printed letter, looked up in
`03_Revised_Panels/PROVENANCE.csv`. The directory is `04_A`.
Script: `create_correlation_heatmap.py`
Marker: this file, **not** `KNOWN_BROKEN.md`. See `RULES.md` rule 1.

## What was superseded

The **size** of the fifty-six row labels. The submitted figure
(`00_GROUND_TRUTH/figures/Figure 4.pdf`) prints them at 3.322 pt; the redraw
prints them at **4.0 pt**. Nothing else moved: same correlation matrix, same
module ordering, same masked diagonal, same colour map, same limits, same
module boundaries, same five module names along the bottom, same fifty-six
names (the five `Mast` rows keep their `Cx_` prefix so that the names stay
distinct).

## History

- 2026-09-11: the author withdrew the labels altogether, on the reasoning
  that 3.3 pt is unreadable and 6 pt would need 130 mm of panel (the rows sit
  at a 4.68 pt pitch on a 92.4 mm axes). The previous version of this file
  recorded that ruling.
- 2026-09-14: the author reversed it. The names go back on the page, at a
  small size, as **the one exemption** to the 6 pt floor in the whole figure
  set. 4.0 pt is the largest type the 4.68 pt pitch holds with clear paper
  between lines. The exemption is declared, with this reason, in
  `12_Figure_Refactor/sweep_panels.py` `FLOOR_EXEMPTIONS`, and
  `sweep_pages.py` holds this panel's region of the page to the same floor.
  Figure 4B, the previous exemption, is at 6 pt since the same date.

## Why there is no Retest-by date

There is nothing to retest. The panel will not reproduce the printed one
again, by design and under the ruling above. The two senses of
`reproduces_published = no` are set out in `RULES.md` rule 1 and enforced by
`verify_panel_provenance.py` checks 3 and 10.

## Record

- Before side: `07_Archive/2026-09-11_round35_published_look_restored/`
- `PROVENANCE.csv`, Figure 4 panel A: `reproduces_published = no`.
