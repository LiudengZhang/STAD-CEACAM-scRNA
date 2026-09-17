# Figure 4 panel E — superseded, not broken

Ruling-date: 2026-09-11
Ruled-by: the author (CIR-26-0753-ET major revision)
Panel: Figure 4 **E** — the printed letter, looked up in
`03_Revised_Panels/PROVENANCE.csv`. The directory is `04_E`, which also holds
panel D.
Script: `create_momac_marker_dotplot.py`
Marker: this file, **not** `KNOWN_BROKEN.md`. See `RULES.md` rule 1.
Applies-to: E

This directory draws **two** printed panels, D and E, from two different
scripts. Only E is superseded; panel D (`create_momac_umap.py`) still
reproduces the printed panel and is untouched. The `Applies-to:` line above is
what keeps this marker off panel D - without it a directory-scoped marker
convicts every panel in the directory, which is `RULES.md` rule 2 biting from
the other side.

## What was superseded

The **absence of a dot-size key**. The submitted figure
(`00_GROUND_TRUTH/figures/Figure 4.pdf`) carries the colour bar and no size
key; this panel now draws one. Nothing else moved: same genes, same seven
MoMac states, same order, same `standard_scale='var'`, same `cmap='Reds'`,
same dots in the same places at the same sizes and colours.

## Why

The dot **area** on this panel encodes the fraction of cells in each group
expressing the gene. With no key, nothing on the page decodes it.

That was known. The script suppressed scanpy's default key with
`legend(show_size_legend=False)`, on the reasoning that a legend the paper had
never printed would be new content and a redraw is visualisation-only
(`RULES.md` rule 3). `PROVENANCE.csv` recorded the consequence and left it
standing, verbatim: "OBSERVATION, recorded and not acted on: the dot area
encodes the fraction of expressing cells in each group, and with no size key
nothing on the page decodes it."

The author's ruling of 2026-09-11 is that **a key which decodes dots already
drawn adds no data, and so is not new content**. The dots do not move; what
changes is that the reader can now read them.

## What was drawn

Five steps — 20, 40, 60, 80, 100% — laid out by
`00_Config/cnsfig/legend.py`, which spaces the circles on a pitch taken from
their own radii rather than at a constant interval. The sizes come from the
DotPlot's own scale (`dot_min`, `dot_max`, `size_exponent`, `largest_dot`), so
the key is drawn at the same scale as the dots it decodes.

The module raises rather than draws if a dotplot varies its dot area and
nothing decodes it. This panel is the reason that check exists.

## Why there is no Retest-by date

There is nothing to retest. The panel will not reproduce the printed one again,
by design and under the ruling above.

## Record

- Before side: `07_Archive/2026-09-11_round34_figure_legibility/`
- `PROVENANCE.csv`, Figure 4 panel E: `reproduces_published = no`.
