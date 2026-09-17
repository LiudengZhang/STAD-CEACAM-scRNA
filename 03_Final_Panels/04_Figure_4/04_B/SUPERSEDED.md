# Figure 4 panel B — superseded, not broken

Ruling-date: 2026-09-14
Ruled-by: the author (CIR-26-0753-ET major revision)
Panel: Figure 4 **B** — the printed letter, looked up in
`03_Revised_Panels/PROVENANCE.csv`. The directory is `04_B`.
Script: `create_merged_networks.py`
Marker: this file, **not** `KNOWN_BROKEN.md`. See `RULES.md` rule 1.

## What was superseded

The **node names** of the five module networks as printed in the submitted
figure (`00_GROUND_TRUTH/figures/Figure 4.pdf`): `C2_CD8`, `C9_CD4`,
`C1_Mono`, `C3_Mac`, `C0_PB`, `C6_B` and their thirty-two siblings. They are
now **one letter and one digit** — `C2`, `H9`, `O1`, `M3`, `P0`, `B6` — the
letter naming the lineage and the digit being the cluster index the state
name already carried. The letters are owned by
`00_Config/shared/labels.py` `NODE_LETTER`, and the Figure 4 legend defines
them in a sentence `04_Manuscript_R1/verify_manuscript_text.py` checks
against that table. All thirty-eight names drawn are distinct; the script
raises otherwise.

Also superseded: the node **circle diameter**, 3.63 mm on the page, is
4.0 mm, so that a 6 pt two-character name sits inside it with clear paper.

Nothing else moved: same five Jaccard matrices, same similarity threshold,
same node selection, same circular seating, same lineage colours, same edge
weights, same module titles.

## Why

At the figure's 6 pt floor the printed names — up to seven characters —
collided round a ring 21.8 mm wide, and on 2026-09-11 they were carried at
3.55 pt as a measured exemption. On 2026-09-14 the author ruled the other
way: the names are shortened to one letter and one digit so they can be set
at 6 pt, and the circles grow to hold them. The exemption is gone from
`12_Figure_Refactor/sweep_panels.py` `FLOOR_EXEMPTIONS`.

## Why there is no Retest-by date

There is nothing to retest. The panel will not reproduce the printed one
again, by design and under the ruling above. See `RULES.md` rule 1 and
`verify_panel_provenance.py` checks 3 and 10.

## Record

- Before side: `07_Archive/2026-09-11_round35_published_look_restored/`
- `PROVENANCE.csv`, Figure 4 panel B: `reproduces_published = no`.
