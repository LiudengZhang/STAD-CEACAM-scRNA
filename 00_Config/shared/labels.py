"""
Authorised label changes between the two drawings of a panel.

A redraw is a visualisation-only change: no number, gene set, cell selection,
statistic or panel content may move. Two changes to *text* are nevertheless
authorised, because they are typography and not content:

  RENAMES   a string printed in a shortened form so that it fits the slot it is
            drawn in - a Hallmark gene-set name, a cell-type name. The value
            behind the label does not move.
  REMOVALS  a string deleted from a panel because it is stated once elsewhere:
            in a shared legend, in the caption, in a shared axis label. The
            value it described is still drawn.

Both are declared here so that `10_Reproduction/compare_panel_content.py`
reports them by name instead of either failing on them or tolerating them
silently. An undeclared string change is still a failure.

Every mapping is read old -> new and is direction-aware. A RENAMES entry
accepts the old string being replaced by the new one and never the reverse, so
a label that drifts back to its long form is still caught. A REMOVALS entry
maps the deleted string to where the reader now finds it: 'legend', 'caption',
'panel title', 'shared y-axis label', and so on.

A worked example, in the shape the sections below take:

    RENAMES_FIGURE_5 = {
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNFA signalling via NF-kB",
    }
    REMOVALS_FIGURE_5 = {
        "Normalised enrichment score": "shared y-axis label",
    }

Sections are keyed by figure and merged at the bottom of this file, so one
figure's declarations can be edited without touching another's. A key declared
with two different values inside one figure's effective table - its own section
together with the supplementary one - is a contradiction and raises on import.
The same key sent to two different destinations by two different figures is
not: two figures can legitimately print "Pre-R" and "Post-R" where both their
scripts set "R". Such a key has no figure-independent destination, so it is
kept out of the merged table and recorded in CROSS_FIGURE_CONFLICTS; see
_merge(). A declaration that never fires is a stale defect, and the gate warns
about it rather than passing quietly.
"""

# --- Figure 2 ---------------------------------------------------------------

RENAMES_FIGURE_2 = {
    # THE STRINGS THE SHIPPED PAGE PRINTS, PUT BACK.
    #
    # Figure 2's panels reach the page as artwork, and an assembler that is not
    # on disk set some of their type. Where that type is still live in
    # Main_Figures/_patched/Figure_2.pdf it can be read off the page, and it is
    # not always what the panel script sets: the page prints "Pre-R" and
    # "Pre-NR" where five panels set "R" and "NR", "Epi_CEACAM5/6" where panel
    # E sets "CEACAM5/6 / Epithelial", "Pre-R vs Pre-NR" where panel F sets
    # "Pre: R vs NR", and an "Epi_" prefix on the nine cluster labels of panel
    # A. The page is the ground truth and the Results text names the page's
    # labels, so the page's string is the one that is drawn and the script's is
    # the one declared here. Not one value moves: every P value the scripts
    # compute already agrees with the page to the printed decimal.
    #
    # The four response labels. Read off the patched page at 4 pt in panels K
    # and L, at 6 pt in E and N.
    "R": "Pre-R",
    "NR": "Pre-NR",
    "R\n(n=4)": "Pre-R",
    "NR\n(n=4)": "Pre-NR",

    # Panel E's title, read off the patched page at 8 pt.
    "CEACAM5/6\nEpithelial": "Epi_CEACAM5/6",

    # Panel F's title, read off the patched page at 8 pt.
    "Pre: R vs NR": "Pre-R vs Pre-NR",

    # NOT DECLARED, AND THE REASON IS MEASURED: the nine cluster labels of
    # panel A. The page prints them at 5.47 pt with an "Epi_" prefix -
    # "Epi_PTMA", "Epi_CEACAM5/6", "Epi_Chief_Like" - and they do not fit that
    # panel at this figure's 6 pt floor. sc.pl.umap locks the axes to equal
    # aspect, so the millimetres between two labels are set by the panel's
    # printed height and no margin recovers any: "Epi_CEACAM5/6" sets 16.23 mm
    # and "Epi_Chief_Like" 14.23 mm, their centres stand 12.87 mm apart, and
    # half their widths sum to 15.23 mm, so the pair overlaps by 2.37 mm. The
    # panel therefore keeps the script's own shorter names, which clear by
    # 1.88 mm, and the difference from the page is reported rather than
    # absorbed. See the header of 02_A/create_epithelial_umap.py.

    # Panel D's annotation. The panel reports two correlations of the same two
    # genes at two levels of pooling, and the figure legend now states the
    # per-cell one and its cell count in words, so the panel prints the
    # metacell value alone rather than the same two numbers twice. Both are
    # still computed and both are still written to
    # 02_D/ceacam_correlation_levels.csv. Declared as a rename rather than a
    # removal because the annotation is one string and only part of it goes:
    # the gate fires a removal only when a whole element is emptied.
    "ρ = 0.72, metacells of 10 (n = 2,168)\nρ = 0.44 per cell (n = 60,937)":
        "ρ = 0.72",

    # The titles of panel L's two boxes. Each box is 22.1 mm wide and about
    # 10 mm of that is the rotated y label and the y tick labels, so the
    # plotting box a title is centred on is about 12 mm; the one-line form sets
    # 18.5 mm at 7 pt and cannot be fitted at all, and on two lines it takes so
    # much of a 20 mm box that the bracket's P value is pushed into it. The
    # title carries the gene alone, as panel K's two titles do and as the page
    # prints panel K's. The compartment is named in the caption, and this
    # panel's y axis already says the values are log2 deconvolved expression.
    "$\\it{CEACAM6}$ (Epi)": "$\\it{CEACAM6}$",
    "$\\it{CEACAM5}$ (Epi)": "$\\it{CEACAM5}$",

    # Panel L's y axis label, re-wrapped. Every word is kept; only the line
    # break is new. A one-line rotated label sets 19.2 mm on a 20.2 mm box, so
    # wherever it is centred it reaches into the 3.8 mm corner the panel letter
    # is printed in. On two lines it sets 12.2 mm and clears the corner.
    "Expression (log2)": "Expression\n(log2)",
}

REMOVALS_FIGURE_2 = {}

# --- Figure 3 ---------------------------------------------------------------

RENAMES_FIGURE_3 = {
    # The two group tick labels of the three paired spatial boxplots H, I and
    # J. At 6 pt "CEACAM-" sets 9.76 mm while the two ticks stand 5.55 mm apart
    # on H, 8.29 on I and 7.29 on J, so the pair overprints on all three; and
    # the labels are themselves what holds the axes narrow, because the fit
    # keeps them inside the canvas, so widening the axes moves them further
    # apart and further out in equal measure. The widest label these panels can
    # carry is about 6.6 mm. The qualifier is the same on every box of all
    # three panels and is stated in the caption; the two groups keep their own
    # names. Same two groups, same order, same tick positions.
    "CEACAM-\nlow": "Low",
    "CEACAM-\nhigh": "High",

    # The titles of the two distance MAPS, K and L. Added 2026-09-10 with the
    # micrometre conversion of those two panels. The published page printed
    # their colour bars as bare numbers - K as 0/2000/4000, L as
    # 500/1000/1500/2000 - in Visium array units, with no unit anywhere on the
    # panel. That was legible while I and J beside them printed
    # "Distance (a.u.)", and stopped being legible when I and J were converted
    # to micrometres earlier the same day: an array unit is 0.2874 um, so the
    # neighbouring panels' unit read across to these two is out by 3.476x.
    #
    # The values are converted by 00_Config/spatial_scale.py's cohort factor
    # and the unit is added to the title rather than to the bar, because
    # REMOVALS_FIGURE_3's stated reason for the bar carrying no label is that
    # the quantity "is named by the panel title above the map". The unit is
    # part of that name, so no new text object is drawn.
    #
    # On a SECOND LINE, measured rather than preferred. "Distance to Stroma
    # (µm)" on one line sets 28.5 mm in a 30.7 mm panel, so its first glyph
    # falls inside the panel-letter keep-out cell and the fit pushes the panel
    # down 2.2 mm to clear it, costing the colour bar a labelled tick. Wrapped,
    # the first line keeps the published title's width and x position to within
    # 0.02 mm, the bar stays at x 25.40 mm where it was published, and it
    # carries four ticks instead of the published three. Same words, one line
    # break, plus the unit.
    "Distance to Stroma": "Distance to Stroma\n(µm)",
    "Distance to Immune": "Distance to Immune\n(µm)",

    # Panel N's y axis label and its title, re-wrapped. Both carry every word
    # they carried before; only where the lines break changes. The axis label
    # is rotated, so its length is vertical: on two lines it sets 35.0 mm
    # against a 36.0 mm panel and cannot clear the panel-letter corner at the
    # top of that column, and on four it sets 17.4 mm. The title sets 45.9 mm
    # on one line while the plotting box is 37 mm wide, and a title is centred
    # on that box, so a title wider than it cannot be fitted at all; on two
    # lines it sets 22.9 mm.
    "Difference in Proportion\n(CEACAM-high - CEACAM-low)":
        "Difference in\nProportion\n(CEACAM-high\n- CEACAM-low)",
    "Immune Recruitment by CEACAM Region":
        "Immune Recruitment\nby CEACAM Region",

    # Two of the eight CD8 state labels on the embedding, panel C. Each label
    # is placed at its state's median embedding position, so two neighbouring
    # states are labelled a few millimetres apart whatever the panel is: at
    # 6 pt "Cytotoxic_CCL" sets 14.5 mm and prints over the label beside it.
    # The states are named by their marker gene, which is how the printed panel
    # names them. Same eight states, same colours, same positions.
    "Cytotoxic_CCL": "CCL",
    "Cytotoxic_DUSP1": "DUSP1",

    # The titles of the four scatter panels A, B, D and E. Each panel draws two
    # axes and each axes titled itself with the cohort it plots, so one cohort
    # name is printed four times over on A and D and four times over on B and
    # E. The cohort, and the sample count where the title carries one, are
    # named once in the caption; each axes keeps its own correlation and its
    # own P value, which are what differ between them. On one line the pair
    # sets 23.6 mm against a plotting box of about 10 mm, and the left axes'
    # title would print into the right axes', so the two statistics take a line
    # each - the same form panels K and M of Figure 5 take.
    #
    # The cohort the caption must name is the page's, not the script's: the
    # shipped page prints "In house (scRNA-seq)" on A and D where the script
    # sets "Primary Cohort (scRNA-seq)", and "TCGA-STAD" on B and E where the
    # script sets "TCGA-STAD (n = 220)" and its three siblings. The page prints
    # no sample count on B or E.
    #
    # A very small P is written out rather than set as a mathtext power of ten:
    # mathtext draws a superscript at 70% of its base size, so an exponent on a
    # 7 pt title prints at 4.9 pt, below the floor this figure set is set to.
    "Primary Cohort (scRNA-seq)\nρ = 0.36, P = 0.04": "ρ = 0.36\nP = 0.04",
    "Primary Cohort (scRNA-seq)\nρ = 0.26, P = 0.1": "ρ = 0.26\nP = 0.1",
    "Primary Cohort (scRNA-seq)\nρ = 0.32, P = 0.07": "ρ = 0.32\nP = 0.07",
    "Primary Cohort (scRNA-seq)\nρ = 0.31, P = 0.08": "ρ = 0.31\nP = 0.08",
    "TCGA-STAD (n = 220)\nρ = 0.32, P = 8$\\times 10^{-7}$":
        "ρ = 0.32\nP = 8e-7",
    "TCGA-STAD (n = 223)\nρ = 0.31, P = 2$\\times 10^{-6}$":
        "ρ = 0.31\nP = 2e-6",
    "TCGA-STAD (n = 384)\nρ = 0.41, P = 9$\\times 10^{-17}$":
        "ρ = 0.41\nP = 9e-17",
    "TCGA-STAD (n = 387)\nρ = 0.42, P = 1$\\times 10^{-17}$":
        "ρ = 0.42\nP = 1e-17",

    # The axis labels of the four scatter panels. Each is about 22 mm tall and
    # carries two square plotting boxes of 7 to 10 mm, with a rotated y label
    # and a column of tick labels in front of each, so the width an axis label
    # may take is the width of its own box and very little more.
    #
    # The x labels carry the gene alone. At 7 pt the unit sets 18.0 mm on
    # panels A and D and 11.0 mm on B and E, against boxes 7 to 10 mm wide
    # standing 14 to 20 mm apart, so the two axes' units print over each other
    # whatever the margins are. The unit is the same on both axes of all four
    # panels and is stated in the caption.
    #
    # The y labels are re-wrapped where the space allows it and shortened where
    # it does not. A rotated label's length is vertical: panel A's sets 35.8 mm
    # on one line against a 23.1 mm panel and 18.0 mm on two; panel D's sets
    # 24.5 mm against 22.9 mm, and two lines would cost the width its shared
    # legend needs, so it keeps the cell subset and the per-cent sign and drops
    # the word the sign carries.
    #
    # The mathtext is written out where it survives. Mathtext draws a subscript
    # or a superscript at 70% of its base size, so either would print at 4.9 pt
    # on a 7 pt label, below the floor this figure set is set to.
    "CEACAM5 (mean log expr.)": "CEACAM5",
    "CEACAM6 (mean log expr.)": "CEACAM6",
    "PD-L1 (CD274) (mean log expr.)": "PD-L1 (CD274)\n(mean log expr.)",
    # Panel D's y label is also a string the script does not print as the page
    # does: the page reads "Tex fraction in CD8+ (%)" and the script sets the
    # same words in a different order. The page's order is restored and one
    # word goes, because the page's own form sets 27.1 mm against a 22.9 mm
    # panel and this label is rotated; on two lines it fits that height but
    # costs 5 mm of width, and at every gutter wide enough to take them the
    # right axes' title prints into the shared legend instead.
    "CD8$^+$ Tex fraction (%)": "Tex in CD8+ (%)",
    "Epi. CEACAM5 (log$_2$+1)": "Epi. CEACAM5",
    "Epi. CEACAM6 (log$_2$+1)": "Epi. CEACAM6",
    "Epi. PD-L1 (CD274) (log$_2$+1)": "Epi. PD-L1\n(CD274)",
    "Tex ssGSEA (z-score)": "Tex ssGSEA\n(z-score)",
}

REMOVALS_FIGURE_3 = {
    # The note naming the test, on the same three panels. At 6 pt
    # "n = 10 samples / Wilcoxon signed-rank (two-sided)" sets 31.6 mm and the
    # widest of the three panels is 29.1 mm, so it cannot be drawn on any of
    # them at the figure's type size. The three panels report one test between
    # them, so it is named once in the caption. The exact two-sided P value
    # stays on each panel, where the reader needs it.
    "n = 10 samples\nWilcoxon signed-rank (two-sided)": "caption",

    # The coordinate furniture of the five spatial maps, F, G, K, L and M. The
    # shipped page prints none of it: each of those panels carries its title,
    # its colour bar and the bar's ticks, and nothing else. Read at 600 dpi to
    # be sure the strings are not drawn as outlines with no text layer - they
    # are not there at all. The coordinates are Visium pixel positions within
    # one section and carry nothing the reader reads off them, and the quantity
    # each bar measures is named by the panel title above the map. The maps,
    # their colours, their scales and their colour-bar ticks are unchanged.
    "X coordinate": "not printed on the published page",
    "Y coordinate": "not printed on the published page",
    "CEACAM Ratio": "the panel title above the map",
    "Epi Density": "the panel title above the map",
    "Dist to Stroma": "the panel title above the map",
    "Dist to Immune": "the panel title above the map",
    "MoMac Score": "the panel title above the map",
}

# --- Figure 4 ---------------------------------------------------------------

# The forest panel's fold-change axis is logarithmic and labels its decades in
# scientific notation. Mathtext draws a superscript at 70% of its base, so the
# exponent of a 6 pt tick label prints at 4.2 pt, below the floor this figure
# is set to. The decades are written out instead. Same tick positions, same
# scale, same limits, same values; only the notation changes. The three
# decades this figure shares with the supplementary cytokine axis are declared
# in the supplementary section and are not repeated here.
RENAMES_FIGURE_4 = {
    "$\\mathdefault{10^{-3}}$": "0.001",
    "$\\mathdefault{10^{1}}$": "10",
    "$\\mathdefault{10^{2}}$": "100",
    "$\\mathdefault{10^{3}}$": "1000",

    # Figure 4 ships as the submitted page plus the label corrections applied
    # to it, so the page prints strings its panel scripts never produced. A
    # drawing that reaches the page without going through that correction has
    # to carry them itself. These are not shortenings: each is a string the
    # printed page carries and the script did not, restored so that the redrawn
    # page prints what the published one prints. No value, row, order or
    # statistic moves with any of them.
    #
    # The C3 monocyte-macrophage state, named on the page in the form the
    # lineage analysis settled on. Restored in panels A, D, E and G.
    "Mac_IL1B": "MoMac_IL1B",
    # The same state naming the axis of both external-cohort boxes, panel H.
    # Set as two lines: at the figure's body size the restored title is 25.9 mm
    # of rotated type against a 23.2 mm box and does not fit on one.
    "Mac_3 Proportion": "MoMac_IL1B\nProportion",
    # The HGNC symbol for the gene the dotplot displayed as CD16, panel E.
    "CD16": "FCGR3A",
    # The two response groups of the post-treatment boxplot, panel C, restored
    # to what the page names them. Another figure declares the same two strings
    # with a different destination, because a tick reading "R" is "Pre-R" on
    # its page and "Post-R" on this one; both are read off the page each is
    # drawn on. The merge records such a key as contested and leaves it out of
    # the figure-independent table, so each figure gets its own answer through
    # BY_FIGURE and neither is given the other's.
    "R": "Post-R",
    "NR": "Post-NR",
    # The panel title the page prints above that boxplot, which the script set
    # not at all. Declared from the empty string because that is what the
    # script drew in the slot the page fills.
    "": "IM-MoMac",
    # The fold-change axis title, panel G, which the page states with the
    # treatment phase of both groups.
    "Fold Change (R/NR)": "Fold Change (Post-R/Post-NR)",
}

REMOVALS_FIGURE_4 = {}

# --- Figure 5 ---------------------------------------------------------------

# The four treatment/response group names on the two regulon panels, D and E.
# One line at 45 degrees needs about twice the width of two horizontal lines,
# and neither panel is 33 mm wide. Same four groups, same order, same tick
# positions; only the line break and the rotation change.
RENAMES_FIGURE_5 = {
    # The longest of the nine Hallmark names of panel A, set with the Greek
    # letters. The gene-set table spells the two symbols out as "alpha" and
    # "k", and the panel printed them that way while every other NF-κB on
    # the page - the four dotplots of panel N, the two axis titles of K and M -
    # printed the Greek. One page named one pathway two ways. The glyphs are
    # drawn from the same embedded family as the rest of the type and are set
    # as literal characters rather than as mathtext, which would draw them at
    # 70% of the base size and put them under the figure's 6 pt floor. Same
    # gene set, same NES, same bar, same order.
    "TNF-alpha Signaling via NF-kB": "TNF-\u03b1 Signaling via NF-\u03baB",

    "Pre-NR": "Pre\nNR",
    "Pre-R": "Pre\nR",
    "Post-NR": "Post\nNR",
    "Post-R": "Post\nR",

    # The thirteen cell types of the NF-kB radar, panel H, and of the cytokine
    # dotplot, panel F. Thirteen names set around a 33 mm circle collide at
    # full length, and the dotplot's label column is otherwise wider than the
    # whole 37 mm panel. Both panels print one vocabulary, so the figure names
    # a cell type one way. Rows, spokes, their order and their values do not
    # move; only the printed form of the name does.
    "Endothelial cells": "Endo",
    "Endothelial": "Endo",
    "Epithelial cells": "Epi",
    "Epithelial": "Epi",
    "Fibroblast": "Fibro",
    "Neutrophils": "Neut",
    "Pericyte": "Peri",
    "Monocytes/Macrophages": "MoMac",
    "Mast cells": "Mast",
    "DC cells": "DC",
    "NK cells": "NK",
    "Plasma cells": "Plasma",
    "CD4+ T cells": "CD4+ T",
    "CD8+ T cells": "CD8+ T",

    # The radar's two series. At full length the legend is more than half the
    # width of the panel and prints over the lower-right spoke names. The two
    # timepoints are named as the group labels of panels D and E name them.
    "Pre-treatment": "Pre",
    "Post-treatment": "Post",

    # The four boxes of panel J. Each is between 18 and 21 mm wide, and the
    # gene is named once in the shared axis label, so a box title carries the
    # cell type alone, in the same short forms the rest of the figure uses.
    "PD-L1 (CD274)\nMonocytes/Macrophages": "MoMac",
    "PD-L1 (CD274)\nEpithelial": "Epi",
    "PD-L1 (CD274)\nFibroblast": "Fibro",
    "PD-L1 (CD274)\nDC cells": "DC",

    # The two correlation panels, K and M. A title is centred on the plotting
    # box, so a title wider than that box cannot be fitted at all - narrowing
    # the box by a millimetre moves the title's edge by half of one. The
    # correlation and its P value therefore take a line each. Same statistics,
    # same values.
    #
    # The exponent is written out rather than set as a mathtext power of ten:
    # mathtext draws a superscript at 70% of its base size, so an exponent on
    # a 7 pt title prints at 4.9 pt, below the floor this figure set is set to.
    "Tex (CD8+)\nρ = 0.74, P = 1$\\times 10^{-6}$":
        "Tex (CD8+)\nρ = 0.74\nP = 1e-6",
    "Th17 (CD4+)\nρ = 0.40, P = 0.02":
        "Th17 (CD4+)\nρ = 0.40\nP = 0.02",

    # A rotated axis label two lines deep costs twice its leading in width, and
    # panels K and L are 22.5 and 18.5 mm wide. The same words, on one line.
    "Tex (CD8+)\nScore": "Tex (CD8+) Score",
    "IL-6/JAK/STAT3\nScore": "IL-6/JAK/STAT3 Score",

    # Panel L names the gene set in full in its axis label, so its title
    # carries the cell type alone, as the box titles of panel J do.
    "IL-6/JAK/STAT3\nCD4+ T cells": "CD4+ T",

    # The response tick labels of the four boxes of panel J and of panel L,
    # with the sample counts taken out and stated in the figure legend. At
    # 6 pt "(n=5)" sets 5.0 mm and the two ticks of an 18 mm box are 2.7 mm
    # apart, so the counts overprint each other wherever they are drawn; even
    # deleting the shared axis label leaves them 0.4 mm short.
    #
    # Declared as renames rather than removals because a tick label is one
    # string and only part of it goes: the gate fires a removal only when a
    # whole element is emptied or absent. The destination is the legend, the
    # same place Figure 2 D's per-cell correlation was moved to.
    "R\n(n=5)": "R",
    "NR\n(n=6)": "NR",
    "NR\n(n=5)": "NR",

    # The eight Hallmark sets of panel N, in their conventional short forms -
    # the same forms the pathway UMAPs of this figure already print, so one
    # figure names a gene set one way. Each is one line: a label set at an
    # angle puts its second line across its neighbour's first, so a two-line
    # name collides whatever its length. The four dotplots each keep their own
    # eight labels; they sort by NES within cell type and print four different
    # orders, so one shared set would name the wrong dot in three of them.
    "TNF-α/NF-κB": "TNFα/NF-κB",
    "Inflammatory\nResponse": "Inflammation",
    "IL-6/JAK/\nSTAT3": "IL-6/STAT3",
    "IFN-γ\nResponse": "IFN-γ",
    "Angio-\ngenesis": "Angiogenesis",

    # Panel N's legend strip is 12.0 mm wide. "Significance" sets 13.2 mm at
    # 7 pt and "Enrichment" 12.4 mm, so neither heading fits its own strip.
    # The three dots carry the thresholds themselves and the colorbar carries
    # the scale, so each heading names the quantity instead.
    "Significance": "P value",
    "Enrichment": "NES",

    # Panel F's two legend headings. scanpy centres them over its legend
    # column, which is 16 mm of this 38.5 mm panel, so at full length they
    # reach back across the dot matrix and are painted over by it: "Fraction of
    # cells" sets 17.2 mm at 7 pt and "Mean expression" 18.8 mm. Same two
    # quantities, named to fit the column they are set in.
    "Fraction of cells\nin group (%)": "Cells in\ngroup (%)",
    "Mean expression\nin group": "Mean expr.\nin group",
}

REMOVALS_FIGURE_5 = {
    # Panel I. The panel script drew four axes - an enrichment curve and a
    # gene-hit strip for each of two cell types - and the published panel
    # prints two. Measured off the submitted page at 600 dpi, the ink inside
    # panel I's rectangle falls in four bands: 81.5-92.2 mm the fibroblast
    # curve, 92.3-94.1 its rank labels, 96.9-107.9 the epithelial curve,
    # 108.1-109.9 its rank labels, with 94.1-96.9 mm empty; and the extracted
    # text holds no "Hits" and no "Gene Rank" anywhere in the rectangle. The
    # published figure is the ground truth, so the two strips are not drawn.
    # Each strip is named here by its own axis label.
    "Hits": "not printed - the published panel draws no gene-hit strip",

    # The two curves are read against one scale, so it is named once, beside
    # both of them. Both curves and every value on them are still drawn, and
    # the two x scales stay separate: the ranked lists are different lengths.
    "Enrichment Score": "shared y-axis label of the two curves",

    # Panel J prints four boxes side by side under one letter. The quantity is
    # named once, beside the leftmost box. The y tick labels stay on all four:
    # the four ranges differ, so one shared limit would redraw the data.
    "PD-L1 (CD274)\nExpression": "shared y-axis label of the leftmost box",

    # Panel N's legend strip, 12.0 mm wide. The colorbar sits directly under
    # the heading that names its scale, so it is not labelled a second time;
    # and the note on the sign of the score is a sentence about how to read the
    # colour, not a value, so it is stated in the caption.
    "NES": "the legend heading above the colorbar",
    "Positive = NR enriched": "caption",
}

# --- Supplementary figures --------------------------------------------------

# The symlog axis of the epithelial-versus-myeloid cytokine panel labels its
# decades in scientific notation. Mathtext draws a superscript at 70% of its
# base size, so an exponent on a 6 pt tick label prints at 4.2 pt - below the
# floor the figure set is set to. The three decades are 0.01, 0.1 and 1, which
# are shorter written out than as powers of ten and put every glyph on that
# axis at the tick size. Same tick positions, same scale, same numbers; only
# the notation is shortened to fit.
RENAMES_SUPPLEMENTARY = {
    "$\\mathdefault{0}$": "0",
    "$\\mathdefault{10^{-2}}$": "0.01",
    "$\\mathdefault{10^{-1}}$": "0.1",
    "$\\mathdefault{10^{0}}$": "1",
}

REMOVALS_SUPPLEMENTARY = {}


#: Keys two figures send to two different destinations, as
#: {kind: {key: {section: destination}}}. See _merge().
CROSS_FIGURE_CONFLICTS = {}


def _merge(kind, sections):
    """One mapping from the per-figure sections.

    A key declared twice with two different values inside ONE figure's
    effective table - its own section together with the supplementary one,
    which is the pair `compare_panel_content.load_aliases` puts in force for a
    run - is a contradiction and raises: that figure would have two authorised
    destinations for one label and no way to choose.

    The same key sent to two different destinations by two DIFFERENT figures is
    not a contradiction, and refusing it was wrong. Figure 2's panels print
    "Pre-R" where their scripts set "R"; Figure 4's print "Post-R" where theirs
    set the same "R". Both are the string the shipped page carries, both are
    read off that page, and each is right for its own figure. The gate is
    always narrowed with `--figure`, and reads `BY_FIGURE` rather than this
    table, so each figure gets its own answer.

    What such a key cannot have is a figure-INDEPENDENT destination. So it is
    left out of the merged table entirely and recorded in
    CROSS_FIGURE_CONFLICTS. A run that omits `--figure` then treats it as
    undeclared and fails on it, which is the conservative outcome; it can never
    be silently given another figure's answer.
    """
    per_figure, merged, origin, conflicts = {}, {}, {}, {}
    supp = dict(sections[-1][1]) if sections and sections[-1][0].endswith(
        "SUPPLEMENTARY") else {}
    for name, mapping in sections:
        for key, value in mapping.items():
            if key in supp and supp[key] != value and not name.endswith(
                    "SUPPLEMENTARY"):
                raise ValueError(
                    f"{kind}: {key!r} is declared as {value!r} in {name} and as "
                    f"{supp[key]!r} in the supplementary section. Those two are "
                    f"in force together, so the label has two authorised "
                    f"destinations and no way to choose; reconcile them.")
            if key in merged and merged[key] != value:
                conflicts.setdefault(key, {})[origin[key]] = merged[key]
                conflicts[key][name] = value
            merged[key] = value
            origin[key] = name
            per_figure.setdefault(key, {})[name] = value
    for key in conflicts:
        merged.pop(key, None)
    if conflicts:
        CROSS_FIGURE_CONFLICTS[kind] = conflicts
    return merged


RENAMES = _merge("RENAMES", [
    ("RENAMES_FIGURE_2", RENAMES_FIGURE_2),
    ("RENAMES_FIGURE_3", RENAMES_FIGURE_3),
    ("RENAMES_FIGURE_4", RENAMES_FIGURE_4),
    ("RENAMES_FIGURE_5", RENAMES_FIGURE_5),
    ("RENAMES_SUPPLEMENTARY", RENAMES_SUPPLEMENTARY),
])

REMOVALS = _merge("REMOVALS", [
    ("REMOVALS_FIGURE_2", REMOVALS_FIGURE_2),
    ("REMOVALS_FIGURE_3", REMOVALS_FIGURE_3),
    ("REMOVALS_FIGURE_4", REMOVALS_FIGURE_4),
    ("REMOVALS_FIGURE_5", REMOVALS_FIGURE_5),
    ("REMOVALS_SUPPLEMENTARY", REMOVALS_SUPPLEMENTARY),
])


#: The same sections, kept apart, so a gate can ask for one figure's
#: declarations. The supplementary section is merged in wherever it is used:
#: supplementary panels are gated by the same tool as the main figures.
BY_FIGURE = {
    "2": (RENAMES_FIGURE_2, REMOVALS_FIGURE_2),
    "3": (RENAMES_FIGURE_3, REMOVALS_FIGURE_3),
    "4": (RENAMES_FIGURE_4, REMOVALS_FIGURE_4),
    "5": (RENAMES_FIGURE_5, REMOVALS_FIGURE_5),
    "S": (RENAMES_SUPPLEMENTARY, REMOVALS_SUPPLEMENTARY),
}
