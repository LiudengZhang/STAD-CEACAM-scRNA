# Figure 2 panels H and I — superseded, not broken

Ruling-date: 2026-09-15
Ruled-by: the author (CIR-26-0753-ET major revision)
Panel: Figure 2 **H** and **I** — the printed letters, looked up in
`03_Revised_Panels/PROVENANCE.csv`. The directory is `02_G`; do not read the
directory as the letter.
Applies-to: H I
Script: `create_mp45_horizontal_boxplot.py`
Marker: this file, **not** `KNOWN_BROKEN.md`. See `RULES.md` rule 1.

## What was superseded

The **P values**. The printed panels annotate the one-tailed exact permutation
test (H: `*`; I: P = 0.08). Under reviewer 1's point R1.3c the test is reported
two-sided everywhere - Results, Table S6, the response letter - and the panels
print those values (P = 0.08 and P = 0.19). The star cannot return: two-sided,
neither program is below 0.05.

## What is back as printed

The **brackets**, since 2026-09-15 (the author's fourth reading). The P bracket
runs from Pre-R to x = 2, the centre of the "ns" bracket over Post-R, Pre-NR
and Post-NR (Kruskal-Wallis homogeneity among the three); read together they
say "one group against these three", which is the design. History: on
2026-09-14 the P bracket was moved to span all four groups (x 0 to 3) and the
"ns" bracket dropped, because a bracket ending at Pre-NR had been read as a
pairwise Pre-R-versus-Pre-NR test (whose P is 0.49 / 0.69, not 0.08 / 0.19);
on 2026-09-15 the "ns" bracket was drawn again and, that evening, the P bracket
was returned to the centre of the three. The Figure 2 legend names both tests.

## What did not move

The scores, the groups, both tests and both P values, at any date.
`compare_panel_content.py` against the round-38 script reports the P bracket's
two right-hand vertices (x 3 -> 2) and the P string's x (1.5 -> 1.0), and
nothing else.
