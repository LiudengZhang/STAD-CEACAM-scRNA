# How a main-figure panel is restyled for Version B

Stage 3 of the cnsplots restyle. This is the recipe every panel under
`11_Version_B_Restyle/Main_Figures/0X_Figure_X/` follows. It is the stage-2
recipe (`Supplementary_New/_restyled/S8_CEACAM_Metaprogram/S8_F/`) applied to
the main figures, with one addition — the mark-rescale factor is now derived
rather than measured, see below.

## The rule

A restyled panel is a **copy** of its Version A script with these changes and
no others:

1. the type comes from `00_Config/panel_style_cns.py` (cnsplots), not from
   `shared.figure_config.use_panel_style()` and not from `N * SCALE`;
2. the canvas is the millimetre box the panel prints in, not four times it;
3. everything specified in points that is *not* type — marker areas, marker
   edge widths, line widths, cap sizes, rules — is multiplied by `MARK` (or
   `AREA` for areas), so its size **relative to the type** is unchanged;
4. margins are millimetres of paper, via `style.margins_mm`;
5. the save is `style.save_panel`, which never passes `bbox_inches`;
6. `style.overflow_mm` is checked before the save.

**Nothing else may change.** Not a path, not a filter, not a threshold, not a
gene list, not a statistic, not a string. The drawing code is the same code.

## The mark rescale factor, derived

Version A drew on a canvas `SCALE` times the printed size and the assembler
fitted it into a millimetre box at some fit `f`:

    printed_A_pt = set_A_pt x SCALE x f

Version B draws 1:1, so `printed_B_pt = set_B_pt`. Holding the panel's
proportions means every length must grow by the same ratio the type grew by:

    ratio = tick_pt / (small_A_pt x SCALE x f)

and a length that was set at `L_A` points on the 4x canvas printed at
`L_A x f`, so the length to set now is

    L_B = L_A x f x ratio = L_A x tick_pt / (small_A_pt x SCALE)

**The assembler fit `f` cancels.** So the factor needs no measurement of the
published page at all:

    MARK = style.tick_pt() / (SMALL_PT * SCALE)
    AREA = MARK ** 2

where `SMALL_PT` is the nominal point size the Version A script used for its
*smallest body type* — the number before the `* SCALE`, almost always the
`fontsize=5 * SCALE` on the tick labels. For the usual `SCALE = 4`,
`SMALL_PT = 5`:

    MARK = 7 / 20 = 0.35        AREA = 0.1225

which is the same 0.35 that stage 1 measured its way to for S8F. Every panel
states `SCALE` and `SMALL_PT` as read out of its own Version A source, so the
factor is auditable per panel.

## The panel box

Start from the panel's **printed rect** in `03_Revised_Panels/panel_rects.csv`
— that is the size the panel occupies on the published page, measured off the
page itself. Then grow it until `style.overflow_mm(fig)` returns all zeros and
the panel is not cramped.

Panels **will** grow. A Figure 5 panel printed at 18.5 x 28.6 mm cannot carry a
y axis label, its tick labels and a title at 7/8 pt in 18.5 mm; ~14 mm of that
width is type before any data is drawn. This is expected and is the same effect
stage 2 absorbed in the supplementary layouts. Each panel records both numbers:

    PRINTED_MM  = the published rect, from panel_rects.csv
    PANEL_W_MM / PANEL_H_MM = the Version B box

## Type

Do not set a font size anywhere unless the panel genuinely needs a size
*relative* to the system, and then take it from `style.tick_pt()` or
`style.body_pt()` — never a literal. Delete `fontsize=N * SCALE`; the rcParams
cnsplots installs already carry 8 pt axis labels/titles, 7 pt ticks and legend,
and 8 pt bold panel letters (drawn by the assembler, not the panel).

`fontweight='bold'` on body text is dropped: cnsplots bolds axis *titles* and
panel letters and nothing else, and following the library rather than the
brief's paraphrase is the standard-methods rule (stage 1 finding 2).

## The gate — every panel, no exceptions

    conda run -n Liudeng_Python_310 python 10_Reproduction/compare_panel_content.py \
        03_Revised_Panels/Main_Figures/<A script> \
        11_Version_B_Restyle/Main_Figures/<B script>

Must print `CONTENT IDENTICAL`. A trailing
`[axis tick density differs ...]` note is allowed and is decoration — the
numeric ruler gained or lost gradations because the axis is a physically
different size. **A CONTENT DIFFERENCE is a stop.** Do not "fix" it by editing
the restyled panel until it agrees; find out which of the two moved, and if the
restyle would change what the panel shows, stop on that panel and report it.

The gate's own sensitivity is proved by
`11_Version_B_Restyle/_audit/gate_positive_control.py`, which must be run
before any negative result from the gate is believed.

## Writes

A restyled script writes only into its own directory under
`11_Version_B_Restyle/`. `03_Revised_Panels/Main_Figures/`,
`_patched/`, `Supplementary_New/_assembled/`, `06_Submission_Package/` and
`00_GROUND_TRUTH/` are frozen: Version A ships today. Nothing is deleted
anywhere; archive by copy.

## AnnData

Read expression from `.raw`, never `.X`. Eight input h5ads carry a double
normalisation from 2026-07-30 that left whole cell rows NaN. Version A's
scripts already do this; do not change how they read.

## Never

A patient or specimen identifier in any file, log or commit message.
