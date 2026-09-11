# Figure 3 panel I — superseded, not broken

Ruling-date: 2026-09-10
Ruled-by: the author (CIR-26-0753-ET major revision)
Panel: Figure 3 **I** — the printed letter, looked up in
`03_Revised_Panels/PROVENANCE.csv`. The directory `03_I` happens to agree here;
four Figure 3 panels do not, so the letter was looked up and not inferred.
Script: `create_spatial_boxplot_stroma.py`
Marker: this file, **not** `KNOWN_BROKEN.md`. See `RULES.md` rule 1.

## What was superseded

The **y axis** of the panel as printed in the submitted figure
(`00_GROUND_TRUTH/figures/Figure 3.pdf`): its label `Distance (a.u.)` and its
tick values 0, 500, 1000, 1500, 2000. The redrawn panel prints
`Distance (µm)` and 0, 200, 400, 600, and does not reproduce the
printed axis. Nothing else moved.

## Why

The panel plotted the `distance_to_stroma` column of `spot_data.csv`
directly. That column is a Euclidean distance in the `x`, `y` of the same
table — full-resolution image pixels — and nothing anywhere in that path
converts it. The panel was therefore right to say `a.u.`.

The Results, however, reported the same quantity as micrometres: "mean
difference 371 µm". Those were array units with a micrometre sign on them.
The figure agreed with the data and the sentence agreed with neither.

Author's ruling 2026-09-10: **convert both sides**, so the reader gets a
physical distance and the page and the paragraph agree. The mean difference
reads 107 µm after conversion, from 371 before.

## The factor, and why there is one and not ten

`00_Config/spatial_scale.py` — the single owner under RULES.md rule 5. The
derivation is the one that was already in the tree, in
`02_New_Analyses/05_R1.6_Spatial_Confounders/scripts/spatial_distance_map.py`,
where it put Figure R2's axes into microns for the response letter and went no
further; it moved rather than being copied, and that script now imports it.

    um_per_unit = 100 um Visium pitch / median nearest-neighbour spacing

Measured over the ten GSE251950 sections: spacing 347.0–354.0 units, factor
0.282476–0.288180 um/unit, a 1.99 per cent spread. **One cohort factor,
0.287351576 um/unit**, is used for all ten. That is a decision and it is
recorded rather than taken quietly:

  * the published quantities are cohort-level, so one constant makes this a
    pure linear rescale — the paired Wilcoxon statistic and the exact P value
    printed on this panel are the ones printed before, to the digit;
  * a per-section factor would re-weight the ten pairs against each other,
    which is an edit to the statistic and not to its unit, and would give
    106.93 um where the cohort factor gives 106.54 um for this panel's
    mean difference, a difference of well under one per cent;
  * `spatial_scale.cohort_um_per_unit()` refuses to return one factor at all
    if the sections ever disagree by more than 5 per cent.

## What moved, exactly

    every y coordinate   x 0.287351576, the same factor everywhere: 23 lines,
                         2 point collections, 2 box patches and the bracket
                         annotation's anchor
    y limit              2300 -> 660.909
    y tick labels        matplotlib's locator on the rescaled limit
    y axis label         'Distance (a.u.)' -> 'Distance (µm)'

    every x coordinate   identical to the last digit
    every other string   identical, including the exact two-sided P value
    canvas               the same millimetres, drawn 1:1 into the same
                         measured printed rectangle

That is measured, not asserted: `10_Reproduction/verify_panel_units.py` runs
`compare_panel_content.py`'s own harvester over the archived pre-conversion
script and this one and checks each of those clauses, and carries nine
mutations — a point moved by one micrometre, the P value changed, the title
shortened, a point moved sideways, a whisker left in array units, the label
given the wrong unit, one group scaled by its own factor, the ticks labelled
in the old units, a tick placed outside the axis — every one of which it is
shown refusing in the same run.

## Why there is no Retest-by date

There is nothing to retest. The panel will not reproduce the printed axis
again, by design and under the ruling above. A `Retest-by:` line here would
file an authorised correction as a defect. The two senses of
`reproduces_published = no` are set out in `RULES.md` rule 1 and enforced by
`03_Revised_Panels/verify_panel_provenance.py` checks 3 and 10.

## Record

- Previous version of the panel and its script:
  `07_Archive/2026-09-10_round22_spatial_distance_microns/`
- The table the Results quote, converted by the same factor at its own source:
  `02_New_Analyses/02_R1.3_TwoSided_Stats_Sweep/outputs/twosided_sweep.csv`,
  checked by `10_Reproduction/verify_distance_units.py`
- `PROVENANCE.csv`, Figure 3 panel I: `reproduces_published = no`, with this
  mechanism in its note.
