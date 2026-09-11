# Figure 3 panel K — superseded, not broken

Ruling-date: 2026-09-10
Ruled-by: the author (CIR-26-0753-ET major revision)
Panel: Figure 3 **K** — the printed letter, looked up in
`03_Revised_Panels/PROVENANCE.csv`. The directory `03_K` happens to agree
here; four Figure 3 panels do not, so the letter was looked up and not
inferred.
Script: `create_spatial_dist_stroma.py`
Marker: this file, **not** `KNOWN_BROKEN.md`. See `RULES.md` rule 1.

## What was superseded

The **colour bar** and the **title** of the panel as printed in the submitted
figure (`00_GROUND_TRUTH/figures/Figure 3.pdf`): a bar reading 0 / 2000 / 4000
with no unit anywhere on the panel, under the title `Distance to Stroma`. The redrawn
panel prints 0 / 500 / 1000 / 1500 under `Distance to Stroma` with `(µm)` on a second line, and
does not reproduce the printed bar.

Nothing else moved. Not one spot changed position and not one spot changed
colour.

## Why

The panel scattered the `distance_to_stroma` column of `spot_data.csv` directly. That
column is a Euclidean distance in the `x`, `y` of the same table — full
resolution image pixels, the Visium array unit — and nothing anywhere in that
path converts it. The published bar was therefore *numerically* right, and it
carried no unit at all, which was survivable only for as long as nothing
nearby claimed one.

On 2026-09-10 panels I and J, the two paired boxplots of the same two
distances, were converted to micrometres under the author's ruling and now
print `Distance (µm)`. That is what made this panel wrong to leave: an array
unit is 0.2874 µm, so a reader who carries the unit off the panel next door
reads every number on this bar as **3.4801 times** what it is.

So the same conversion is applied here, by the same factor and from the same
owner.

## The factor

`00_Config/spatial_scale.py` — the single owner under `RULES.md` rule 5,
imported and not re-derived.

    um_per_unit = 100 µm Visium pitch / median nearest-neighbour spacing
                = 0.287351576 µm per array unit   (one cohort factor)

It is taken off the **whole** ten-section spot table and not off this
section's slice, so this map's unit is the same unit panels I and J print.

## Why the unit is in the title, and on a second line

`REMOVALS_FIGURE_3`'s stated reason for this bar carrying no label is that the
quantity "is named by the panel title above the map". The unit is part of that
name, so putting it there adds no new text object to the panel.

The second line is a **measurement, not a preference**. Three layouts were
drawn and the ink read back off the PDFs:

    published                    title x 3.13-25.64 mm, y 1.91-4.38
                                 bar at x 25.40 mm, ticks 0 / 2000 / 4000
    "Distance to Stroma (µm)"      one line   title x 0.60-29.17, y 4.12-6.59
                                 ticks 0 / 1000
    "µm" on the bar              title x 0.60-23.10, y 4.02-6.49
                                 ticks 0 / 1000
    "Distance to Stroma\n(µm)"    ADOPTED    title x 3.11-25.62, y 0.52-2.99
                                 bar at x 25.40 mm, ticks 0 / 500 / 1000 / 1500

A one-line title widens past the panel, so its first glyph lands inside the
panel-letter keep-out cell and `fit_margins` pushes the whole panel down by
2.2 mm to clear it — taking 2.2 mm off the colour bar, which then carries two
labelled ticks instead of three or four. Labelling the bar itself does the
same thing from the other side. The two-line title keeps the first line's
width, and therefore its x position, within 0.02 mm of the published page,
keeps the bar where it was, and leaves it more ticks than the page had.

## What moved, exactly

    every mapped value   x 0.287351576, the same factor everywhere, all 2,419
    colour bar limit     5857.182 -> 1683.070
    colour bar ticks     matplotlib's locator on the new limit
    title                'Distance to Stroma' -> 'Distance to Stroma' + a second line '(µm)'

    every spot's x       identical to the last digit
    every spot's y       identical to the last digit
    every spot's colour  identical — `colour_ranks` equal over all 2,419 spots.
                         `Normalize` maps the data's own range onto the colour
                         map, so one positive constant leaves every normalised
                         position where it was.
    alpha, marker size   identical
    both axis limits     identical
    both rulers          identical (this map prints neither)
    every other string   identical
    canvas               the same millimetres, drawn 1:1 into the same measured
                         printed rectangle 4.6,88.2,35.3,114.6 mm

That is measured, not asserted. `10_Reproduction/verify_panel_units_maps.py`
runs `compare_panel_content.py`'s own harvester over the archived
pre-conversion script and this one, checks each clause above across **both**
axes of the panel — the map and its colour bar — and carries twelve
mutations, every one of which it is shown refusing in the same run:

    a spot moved sideways          the bar left in array units
    a spot moved up                the bar scaled by a second factor
    one value left in array units  the ticks relabelled by hand
    one value nudged by 1 µm       one tick label left as published
    the array scaled by 1.005      the bar widened
    the unit dropped from the title
    two spots' colours swapped

`10_Reproduction/verify_localised_change.py` at 600 dpi shows that all 355,231
pixels that differ on the reassembled page lie inside the printed rectangles
of K and L, and none outside.

## Why there is no Retest-by date

There is nothing to retest. The panel will not reproduce the printed colour
bar again, by design and under the ruling above. A `Retest-by:` line here
would file an authorised correction as a defect. The two senses of
`reproduces_published = no` are set out in `RULES.md` rule 1 and enforced by
`03_Revised_Panels/verify_panel_provenance.py` checks 3 and 10.

## Record

- Previous version of the panel and its script:
  `07_Archive/2026-09-10_round23_fig3KL_distance_microns/`
- The companion conversion of panels I and J, one day earlier in the same
  ruling: `03_I/SUPERSEDED.md` and `03_J/SUPERSEDED.md`
- `PROVENANCE.csv`, Figure 3 panel K: `reproduces_published = no`, with this
  mechanism in its note.

## Still open, and not this panel's to close

Supplementary Figure S4 panels D and E are the all-ten-section versions of
these same two maps, read from the same two columns of the same
`spot_data.csv`, and their colour bars are **also** in array units — labelled
`Distance to stroma` and `Distance to immune`, with no unit. They are in the
same position this panel was in and need the same treatment. They did not get
it here, and the reason is recorded rather than left implicit: the shipped
`S4_Spatial_Validation.pdf` is not built from those scripts. It is carried
over from `01_Reviewer_Materials/figures_submitted/`, and its text layer shows
it was edited by hand in Adobe Illustrator 30.2 after assembly — the five
`*GC6-PM = Peritoneal metastasis (paired with GC6 primary)` footnotes the
panel scripts draw are absent from it. Rebuilding S4 to carry the unit would
also reinstate those footnotes, which is a separate question and one for the
author.
