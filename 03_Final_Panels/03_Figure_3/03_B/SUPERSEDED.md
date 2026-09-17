# Figure 3 panel D — superseded, not broken

Ruling-date: 2026-09-11
Ruled-by: the author (CIR-26-0753-ET major revision)
Panel: Figure 3 **D** — the printed letter, looked up in
`03_Revised_Panels/PROVENANCE.csv`. The directory is `03_B`; do not read the
directory as the letter.
Script: `create_ceacam_tex_scatter.py`
Marker: this file, **not** `KNOWN_BROKEN.md`. See `RULES.md` rule 1.

## What was superseded

The **duplicate y-axis label on the right sub-panel**. The panel holds two
scatters, of the same quantity against CEACAM5 and against CEACAM6, and each
carried its own `Tex in CD8+ (%)`. The label is now set once, on the left. Both
sub-panels keep their own ticks. Nothing else moved: same points, same colours,
same fits, same rho and P values, same axis limits, same legend.

## Why

The two labels collided, and the obvious fix made it worse.

The right sub-panel's rotated y label sits in the gutter between the two
plotting boxes. The left sub-panel's x label, `CEACAM5`, is **wider than the
box it is centred on**, so it overhangs into that same gutter from the other
side. Measured on the shipped panel the two came within **0.13 pt** and printed
as one string.

Widening the gutter is self-defeating: `wspace` is a fraction of the axes
width, so opening it shrinks each plotting box, which pushes the fixed-width x
label further out. Tried and measured — at `wspace=1.05` the two **touched, at
0.00 pt**, worse than the 0.13 pt it started at.

Naming the shared quantity once empties the gutter instead of fighting over it.

## What nothing caught it

`check_restyled_panel.py` tested for **shared area** between two strings, and
two strings that abut share none. It reported `overlaps 0 text pair(s)` and
`PASS` on this panel while it printed as one word. The check now measures
**clearance** and convicts the old drawing at 0.13 pt; it ships nine controls,
including a negative one so the exemption for a wrapped title cannot become a
hole. `RULES.md` rule 4.

## Why there is no Retest-by date

There is nothing to retest. The panel will not reproduce the printed one again,
by design and under the ruling above.

## Record

- Before side: `07_Archive/2026-09-11_round34_figure_legibility/`
- `PROVENANCE.csv`, Figure 3 panel D: `reproduces_published = no`.
