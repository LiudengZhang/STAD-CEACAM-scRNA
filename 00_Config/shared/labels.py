"""
Authorised label changes between the two drawings of a panel.

REVISED 2026-09-14 (and 2026-09-11). Most of the shortenings this file used to declare were
reversed when the panels were given the blank bottom of the page
(12_Figure_Refactor/build_grid_v2.py). Where a published string is back in
full its entry is gone; where it is back but wrapped, the entry now maps to
the wrapped form. The ones that remain are the ones the arithmetic still
forces, and each panel script records why beside its own code.

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

# --- Figure 4 B: the node letters -------------------------------------------
#
# The author's ruling of 2026-09-14: a network node is named by one letter
# and one digit, so the name sits inside its circle at 6 pt. The digit is the
# cluster index the state name carries; the letter is this table's, keyed by
# the lineage token of the state name (the second underscore field). The
# Figure 4 legend defines the letters, and 04_Manuscript_R1/
# verify_manuscript_text.py checks that sentence against this table - one
# owner, two readers. CD4 and CD8 share an initial, as do Mono/Mac/Mast and
# NK/NKT/Neu, which is why the letters are not initials throughout.
NODE_LETTER = {
    "CD4": "H",       # helper T
    "CD8": "C",       # cytotoxic T
    "NKT": "K",
    "NK": "K",
    "Prolif": "K",    # the one proliferating state in IM-T/NK/DC is NK
    "Mono": "O",
    "Mac": "M",
    "MoMac": "M",
    "Mast": "S",
    "DC": "D",
    "Neu": "N",
    "B": "B",
    "Plasma": "P",    # plasmablast
}
#: The definition the Figure 4 legend prints, generated from NODE_LETTER so the
#: two cannot drift. Read by 01_Main_Text/edits.py and verify_manuscript_text.py.
NODE_LETTER_LEGEND = ("H, CD4+ T; C, CD8+ T; K, NK/NKT; O, monocyte; "
                      "M, macrophage; S, mast; D, dendritic; N, neutrophil; "
                      "B, B cell; P, plasmablast")
_LEGEND_LETTERS = {tok.split(",")[0].strip() for tok in NODE_LETTER_LEGEND.split(";")}
if _LEGEND_LETTERS != set(NODE_LETTER.values()):
    raise ValueError(f"NODE_LETTER_LEGEND names {sorted(_LEGEND_LETTERS)} but "
                     f"NODE_LETTER uses {sorted(set(NODE_LETTER.values()))}")

# --- Gene symbols: italic, digits included ----------------------------------
#
# The author's ruling of 2026-09-14 (evening): a gene symbol prints in italic
# wherever it names a gene - the digit too - and nowhere else. This is the
# list of every symbol Figures 2-5 print. 10_Reproduction/
# check_restyled_panel.py:check_italics reads the faces back off the page and
# holds it in both directions: a symbol here set in the regular face fails,
# and an italic string that is not here fails, so the list is exactly what is
# italic on the page. A symbol that is part of a cluster name (Mac_TREM2,
# MoMac_IL1B, Epi_CEACAM5/6 - preceded or followed by an underscore or a
# slash) names a cluster, not a gene, and is upright. The two IHC panels
# (2M, 2N) name proteins and are upright by declaration in the gate.
GENE_SYMBOLS = frozenset({
    # Figure 2
    "CEACAM5", "CEACAM6", "CD274",
    "PDZK1IP1", "PHGR1", "APOL1", "SAA1", "DUOX2", "MUC17", "HSPA1A",
    "EDN1", "SPINK1", "FXYD3", "SEMA3B", "SEC61G", "DNAJB1", "MUC4",
    "IDO1", "CD70", "NT5E", "CD47", "PVRIG", "CEACAM1", "VSIR", "CD276",
    "CD96", "CD27", "ADORA2A", "ENTPD1", "TNFRSF9", "SIGLEC15", "TNFRSF18",
    "HAVCR2", "ICOS", "BTLA", "KLRC1", "CD40", "TIGIT", "TNFRSF4", "LAG3",
    "CTLA4", "PDCD1", "VTCN1",
    # Figure 4
    "TREM2", "C1QA", "APOE", "CD14", "S100A8", "VCAN", "HLA-DRA", "CST3",
    "CLEC10A", "IL1B", "TNF", "CXCL8", "FCGR3A", "CDKN1C", "LST1", "MKI67",
    "TOP2A", "STMN1", "MT1G", "MT2A", "MT1X",
    # Figure 5
    "BACH1", "NFKB1", "IL6", "IL1A", "ADAM17", "TGFB1", "CCL2", "GLG1",
    "TNFSF15", "ITGAM", "IL1RN", "ANXA1", "ICAM1", "PIK3CB", "CXCL5",
})

# --- Figure 1 ---------------------------------------------------------------
# Redrawn at 1:1 into slots on 2026-09-15 (the page was a 254 mm landscape
# whose type printed at two thirds). Panel B's twelve on-plot names are the
# published page's, read off 00_GROUND_TRUTH/figures/Figure 1.pdf; the script
# that had drawn the panel set shorter forms, which never reached the page.
RENAMES_FIGURE_1 = {
    "Epithelial": "Epithelial Cells",
    "T/NK": "T/NK Cells",
    "Mono/Mac": "Monocytes/Macrophages",
    "Plasma": "Plasma Cells",
    "B cells": "B Cells",
    "Endothelial": "Endothelial Cells",
    "Mast": "Mast Cells",
    "DC": "Dendritic Cells",
}
REMOVALS_FIGURE_1 = {}

# --- Figure 2 ---------------------------------------------------------------

RENAMES_FIGURE_2 = {
    # ROUND 35 -> ROUND 36 (2026-09-14): strings the 2026-09-11 drawing had
    # shortened or wrapped and the 2026-09-14 drawing prints in the page's
    # form again. Read old -> new by compare_panel_content.py when the
    # round-35 script is the old side; a Version A script never sets the
    # left-hand form, so none of these fires against it.
    '$\\it{CEACAM6}$': '$\\it{CEACAM6}$ (Epi)',
    '$\\it{CEACAM5}$': '$\\it{CEACAM5}$ (Epi)',

    # P values at two decimals since 2026-09-14 (panel_style_cns.p_label): the
    # page's own form. The value behind each is unchanged and verify_numbers
    # checks it at three decimals against its table.
    "P = 0.057": "P = 0.06",
    "P = 0.083": "P = 0.08",
    "P = 0.188": "P = 0.19",
    "P = 0.114": "P = 0.11",
    "P = 0.056": "P = 0.06",
    "P = 0.059": "P = 0.06",

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
    # "(Epi)" is back in both titles since 2026-09-14; no entry.

    # Panel L's y axis label, re-wrapped. Every word is kept; only the line
    # break is new. A one-line rotated label sets 19.2 mm on a 20.2 mm box, so
    # wherever it is centred it reaches into the 3.8 mm corner the panel letter
    # is printed in. On two lines it sets 12.2 mm and clears the corner.
    "Expression (log2)": "Expression\n(log2)",
}

REMOVALS_FIGURE_2 = {
    # 2026-09-14: strings that moved rather than vanished.
    'Expression': "the left box's y label (one label per pair since 2026-09-14)",
    'Expression\n(log2)': "the left box's y label (one label per pair since 2026-09-14)",
    '$\\it{CEACAM5}$': "drawn as canvas text above the right box of L, as '$\\it{CEACAM5}$ (Epi)' (2026-09-14)",
}

# --- Figure 3 ---------------------------------------------------------------

RENAMES_FIGURE_3 = {
    # 2026-09-14 (evening): the scatter pairs are 13.5 mm squares now, and
    # the one-line y label of D set 21.6 mm; wrapped, every word kept.
    "CD8+ Tex fraction (%)": "CD8+ Tex\nfraction (%)",
    # The dataset titles lose the parenthetical and the statistics line:
    # rho is drawn inside the box and P goes to the legend (author's ruling
    # of the same date; cnsfig.corr_stats is the owner of the P values).
    "In house\n(scRNA-seq)": "In house",
    # ROUND 35 -> ROUND 36 (2026-09-14): strings the 2026-09-11 drawing had
    # shortened or wrapped and the 2026-09-14 drawing prints in the page's
    # form again. Read old -> new by compare_panel_content.py when the
    # round-35 script is the old side; a Version A script never sets the
    # left-hand form, so none of these fires against it.
    'Low': 'CEACAM-\nlow',
    'High': 'CEACAM-\nhigh',
    'In house\n(scRNA-seq)\nρ = 0.36\nP = 0.04': 'In house\n(scRNA-seq)\nρ = 0.36, P = 0.04',
    'In house\n(scRNA-seq)\nρ = 0.26\nP = 0.1': 'In house\n(scRNA-seq)\nρ = 0.26, P = 0.1',
    'In house\n(scRNA-seq)\nρ = 0.32\nP = 0.07': 'In house\n(scRNA-seq)\nρ = 0.32, P = 0.07',
    'In house\n(scRNA-seq)\nρ = 0.31\nP = 0.08': 'In house\n(scRNA-seq)\nρ = 0.31, P = 0.08',
    'TCGA-STAD\n(n = 220)\nρ = 0.32\nP = 8e-7': 'TCGA-STAD\nρ = 0.32, P = 8e-7',
    'TCGA-STAD\n(n = 223)\nρ = 0.31\nP = 2e-6': 'TCGA-STAD\nρ = 0.31, P = 2e-6',
    'TCGA-STAD\n(n = 384)\nρ = 0.41\nP = 9e-17': 'TCGA-STAD\nρ = 0.41, P = 9e-17',
    'TCGA-STAD\n(n = 387)\nρ = 0.42\nP = 1e-17': 'TCGA-STAD\nρ = 0.42, P = 1e-17',

    # P values at two decimals since 2026-09-14 (p_label); values unchanged.
    "P = 0.014": "P = 0.01",
    "P = 0.084": "P = 0.08",
    "P = 0.064": "P = 0.06",

    # The two group tick labels of the three paired spatial boxplots H, I and
    # J. At 6 pt "CEACAM-" sets 9.76 mm while the two ticks stand 5.55 mm apart
    # on H, 8.29 on I and 7.29 on J, so the pair overprints on all three; and
    # the labels are themselves what holds the axes narrow, because the fit
    # keeps them inside the canvas, so widening the axes moves them further
    # apart and further out in equal measure. The widest label these panels can
    # carry is about 6.6 mm. The qualifier is the same on every box of all
    # three panels and is stated in the caption; the two groups keep their own
    # names. Same two groups, same order, same tick positions.
    # The two-line group names are back since 2026-09-14 (the axis is widened
    # at its ends instead); no entry.

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

    # Two of the eight CD8 state labels on the embedding, panel C. Each label
    # is placed at its state's median embedding position, so two neighbouring
    # states are labelled a few millimetres apart whatever the panel is: at
    # 6 pt "Cytotoxic_CCL" sets 14.5 mm and prints over the label beside it.
    # The states are named by their marker gene, which is how the printed panel
    # names them. Same eight states, same colours, same positions.
    "Cytotoxic_CCL": "TCD8_CCL",
    "Cytotoxic_DUSP1": "TCD8_DUSP1",

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
    # Two lines since 2026-09-14, as the page sets them: the dataset on one,
    # the statistics on the other.
    "Primary Cohort (scRNA-seq)\nρ = 0.36, P = 0.04": "In house (scRNA-seq)\nρ = 0.36, P = 0.04",
    "Primary Cohort (scRNA-seq)\nρ = 0.26, P = 0.1": "In house (scRNA-seq)\nρ = 0.26, P = 0.1",
    "Primary Cohort (scRNA-seq)\nρ = 0.32, P = 0.07": "In house (scRNA-seq)\nρ = 0.32, P = 0.07",
    "Primary Cohort (scRNA-seq)\nρ = 0.31, P = 0.08": "In house (scRNA-seq)\nρ = 0.31, P = 0.08",
    "TCGA-STAD (n = 220)\nρ = 0.32, P = 8$\\times 10^{-7}$": "TCGA-STAD\nρ = 0.32, P = 8e-7",
    "TCGA-STAD (n = 223)\nρ = 0.31, P = 2$\\times 10^{-6}$": "TCGA-STAD\nρ = 0.31, P = 2e-6",
    "TCGA-STAD (n = 384)\nρ = 0.41, P = 9$\\times 10^{-17}$": "TCGA-STAD\nρ = 0.41, P = 9e-17",
    "TCGA-STAD (n = 387)\nρ = 0.42, P = 1$\\times 10^{-17}$": "TCGA-STAD\nρ = 0.42, P = 1e-17",

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
    "CD8$^+$ Tex fraction (%)": "CD8+ Tex fraction (%)",
    "Epi. CEACAM5 (log$_2$+1)": "Epi. CEACAM5",
    "Epi. CEACAM6 (log$_2$+1)": "Epi. CEACAM6",
    "Epi. PD-L1 (CD274) (log$_2$+1)": "Epi. PD-L1\n(CD274)",
    "Tex ssGSEA (z-score)": "Tex ssGSEA\n(z-score)",
}

REMOVALS_FIGURE_3 = {
    # 2026-09-14 (evening), the author's ruling: the scatter titles carry the
    # dataset alone; rho is drawn inside each box and P goes to the legend,
    # generated from cnsfig.corr_stats. The eight round-36 titles, by name.
    "In house\n(scRNA-seq)\nρ = 0.36, P = 0.04": "rho inside the box; P in the legend",
    "In house\n(scRNA-seq)\nρ = 0.26, P = 0.1": "rho inside the box; P in the legend",
    "In house\n(scRNA-seq)\nρ = 0.32, P = 0.07": "rho inside the box; P in the legend",
    "In house\n(scRNA-seq)\nρ = 0.31, P = 0.08": "rho inside the box; P in the legend",
    "TCGA-STAD\nρ = 0.32, P = 8e-7": "rho inside the box; P in the legend",
    "TCGA-STAD\nρ = 0.31, P = 2e-6": "rho inside the box; P in the legend",
    "TCGA-STAD\nρ = 0.41, P = 9e-17": "rho inside the box; P in the legend",
    "TCGA-STAD\nρ = 0.42, P = 1e-17": "rho inside the box; P in the legend",
    # 2026-09-14 (evening): A's y label drops its unit line; the legend gives
    # the unit for both axes of A and D.
    "(mean log expr.)": "caption",
    # B and E: the x labels drop "Epi." - on 13.5 mm boxes 1.5 mm apart the
    # two 15 mm labels overprinted; the y label and the legend say the
    # expression is epithelial (BayesPrism-deconvolved).
    "Epi. CEACAM5": "caption",
    "Epi. CEACAM6": "caption",
    # 2026-09-14: the qualifier is back inside the tick labels.
    'CEACAM region': 'the two tick labels, CEACAM-low and CEACAM-high, since 2026-09-14',

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
    # Panel G's x label, on two lines since 2026-09-15 (the author's fourth
    # reading): on one line it is 33 mm of type centred on a 15 mm frame and
    # was holding 9 mm of blank on either side of the forest plot. Same
    # words; a line break.
    "Fold Change (Post-R/Post-NR)": "Fold Change\n(Post-R/Post-NR)",
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
    # 2026-09-14 (evening): panel I's y axes each carry the short form; the
    # legend reads "enrichment score (ES)".
    "Enrichment Score": "ES",
    # ROUND 35 -> ROUND 36 (2026-09-14): strings the 2026-09-11 drawing had
    # shortened or wrapped and the 2026-09-14 drawing prints in the page's
    # form again. Read old -> new by compare_panel_content.py when the
    # round-35 script is the old side; a Version A script never sets the
    # left-hand form, so none of these fires against it.
    'Pre\nNR': 'Pre-NR',
    'Pre\nR': 'Pre-R',
    'Post\nNR': 'Post-NR',
    'Post\nR': 'Post-R',
    'Pre': 'Pre-treatment',
    'Post': 'Post-treatment',
    'MoMac': 'Monocytes/\nMacrophages',
    'CD4+ T': 'IL-6/JAK/STAT3\nCD4+ T cells',
    'R': 'R\n(n=5)',
    'NR': 'NR\n(n=6)',
    # Panel J's leftmost y label, one line at the tick size since 2026-09-14.
    'PD-L1 (CD274)\nExpression': 'PD-L1 (CD274) Expression',

    # Panel A's x label, the page's form restored 2026-09-14 (the script set
    # the short form on 2026-09-11).
    "NES": "NES (Post-NR/Post-R)",

    # P values at two decimals since 2026-09-14 (p_label); 0.052 and 0.055 keep
    # three because two would print the borderline '0.05'. Values unchanged.
    "P = 0.924": "P = 0.92",
    "P = 0.035": "P = 0.04",
    "P = 0.898": "P = 0.90",
    # 2026-09-16 (the author's fifth reading: "no need for the actual P
    # values, ns is enough, and not in grey"): the Kruskal-Wallis bracket
    # over the other three groups on D and E prints "ns" in black, as 2H/I
    # print theirs (cnsfig.boxes.bracket kind="omnibus"). The test and its
    # values (0.924, 0.898) are unchanged; the legend names the test.
    "P = 0.92": "ns",
    "P = 0.90": "ns",
    "P = 0.038": "P = 0.04",
    "P = 0.082": "P = 0.08",
    "P = 0.032": "P = 0.03",
    "P = 0.030": "P = 0.03",

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

    # One line at 45 degrees since 2026-09-14, as the page sets them; no
    # entry.

    # The thirteen cell types of the NF-kB radar, panel H, and of the cytokine
    # dotplot, panel F. Thirteen names set around a 33 mm circle collide at
    # full length, and the dotplot's label column is otherwise wider than the
    # whole 37 mm panel. Both panels print one vocabulary, so the figure names
    # a cell type one way. Rows, spokes, their order and their values do not
    # move; only the printed form of the name does.
    # Since 2026-09-14 (evening) panel F prints all thirteen on one line:
    # the row-name column is 23 mm and the key sits under the matrix.
    # The radar's five short spoke names (Endo, Epi, Fibro, Neut, Peri) are
    # gone since 2026-09-14: the ring is drawn smaller so the full names fit.

    # The radar's two series. At full length the legend is more than half the
    # width of the panel and prints over the lower-right spoke names. The two
    # timepoints are named as the group labels of panels D and E name them.
    # The radar's series key prints the two names in full since 2026-09-14.

    # The four boxes of panel J. Each is between 18 and 21 mm wide, and the
    # gene is named once in the shared axis label, so a box title carries the
    # cell type alone, in the same short forms the rest of the figure uses.
    # Since 2026-09-14 the cell types are named in full, the widest wrapped.
    "PD-L1 (CD274)\nMonocytes/Macrophages": "Monocytes/\nMacrophages",
    "PD-L1 (CD274)\nEpithelial": "Epithelial",
    "PD-L1 (CD274)\nFibroblast": "Fibroblast",
    "PD-L1 (CD274)\nDC cells": "DC cells",

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
    # Panel L's title is the page's two lines again since 2026-09-14.

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
    # The counts are back in the tick labels since 2026-09-14 (the axis is
    # widened at its ends instead); no entry.

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

    # Panel F's two legend headings. scanpy centres them over its legend
    # column, which is 16 mm of this 38.5 mm panel, so at full length they
    # reach back across the dot matrix and are painted over by it: "Fraction of
    # cells" sets 17.2 mm at 7 pt and "Mean expression" 18.8 mm. Same two
    # quantities, named to fit the column they are set in.
    "Fraction of cells\nin group (%)": "Fraction of\ncells in\ngroup (%)",
    "Mean expression\nin group": "Mean\nexpression\nin group",
}

REMOVALS_FIGURE_5 = {
    # 2026-09-14 (evening): K and M, as Figure 3's scatters - rho inside, P in
    # the legend (cnsfig.corr_stats).
    "Tex (CD8+)\nρ = 0.74\nP = 1e-6": "rho inside the box; P in the legend",
    "Th17 (CD4+)\nρ = 0.40\nP = 0.02": "rho inside the box; P in the legend",
    # 2026-09-14: the radar's thirteen spoke names are no longer tick labels.
    # They are drawn as text along their spokes, in full (the five short forms
    # of 2026-09-11 - Endo, Epi, Fibro, Neut, Peri - become the page's words).
    'B cells': 'drawn as text along its spoke, in full (2026-09-14)',
    'DC': 'drawn as text along its spoke, in full (2026-09-14)',
    'Endo': 'drawn as text along its spoke, in full (2026-09-14)',
    'Epi': 'drawn as text along its spoke, in full (2026-09-14)',
    'Fibro': 'drawn as text along its spoke, in full (2026-09-14)',
    'Mast': 'drawn as text along its spoke, in full (2026-09-14)',
    'MoMac': 'drawn as text along its spoke, in full (2026-09-14)',
    'Neut': 'drawn as text along its spoke, in full (2026-09-14)',
    'NK': 'drawn as text along its spoke, in full (2026-09-14)',
    'Peri': 'drawn as text along its spoke, in full (2026-09-14)',
    'Plasma': 'drawn as text along its spoke, in full (2026-09-14)',
    'CD4+ T': 'drawn as text along its spoke, in full (2026-09-14)',
    'CD8+ T': 'drawn as text along its spoke, in full (2026-09-14)',

    # 2026-09-14: titles that moved onto the canvas because their plotting box
    # could not centre them.
    'MoMac': "drawn as canvas text above panel J's first box, 'Monocytes/\\nMacrophages' (2026-09-14)",
    'CD4+ T': "drawn as canvas text above panel L, 'IL-6/JAK/STAT3\\nCD4+ T cells' (2026-09-14)",

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
    # 2026-09-15: the S7 brackets print their P through panel_style_cns.p_label
    # (two decimals, a star below 0.05), as every main-figure bracket does; the
    # analysis modules' own _panel functions set three decimals. Read old ->
    # new by the drivers' --check.
    "P = 0.343": "P = 0.34",
    "P = 0.057": "P = 0.06",
    "P = 0.200": "P = 0.20",
    "P = 0.486": "P = 0.49",
    "P = 0.063": "P = 0.06",
    "P = 1.000": "P = 1.00",
    "P = 0.662": "P = 0.66",
    "P = 0.686": "P = 0.69",
    "P = 0.286": "P = 0.29",
    "P = 0.762": "P = 0.76",
    "P = 0.537": "P = 0.54",
    # 2026-09-16, annotated supplementary proof: the composite is named
    # concisely on the panel; its definition remains in the legend.
    "Summed (published)": "Summed",
    "$\\mathdefault{0}$": "0",
    "$\\mathdefault{10^{-2}}$": "0.01",
    "$\\mathdefault{10^{-1}}$": "0.1",
    "$\\mathdefault{10^{0}}$": "1",
}

REMOVALS_SUPPLEMENTARY = {
    # 2026-09-16, annotated supplementary proof: sample size and test belong
    # in the figure legend rather than as a second caption above panel S8C.
    "Immunohistochemistry, n = 4 responders vs 4 non-responders (two-sided Mann-Whitney)":
        "the Figure S8 legend states the sample size and test",
}


# --- Supplementary Figures S1-S6, redrawn at 1:1 on 2026-09-15 (evening) ---
# The six submitted pages were carried over untouched until the author's
# fourth reading ("the type in the supplementary figures still looks off");
# every panel is now redrawn into Supplementary_New/S<n>_*/ from a copy of
# its predecessor script. These are the string changes those redraws make,
# read old -> new by compare_panel_content.py. They are in force for the
# supplementary gate ONLY (BY_FIGURE["S"]), not merged into the main
# figures' tables: S2 E and Figure 2 E draw the same predecessor title to two
# different printed titles, which is a cross-figure difference, not a
# contradiction (see _merge).
RENAMES_S1_S6 = {
    # S2 D and E stand side by side and each names its tissue, as the printed
    # page set them (00_GROUND_TRUTH/figures/S2_*.pdf: "Primary Tumor,
    # Epi_CEACAM5/6 / (Tumor-score adjusted)" and "Liver Met.,
    # Epi_CEACAM5/6"); the two predecessor scripts set titles that did not.
    "CEACAM5/6 Epithelial\n(tumor-score adjusted)": "Primary Tumor, Epi_CEACAM5/6\n(Tumor-score adjusted)",
    "CEACAM5/6\nEpithelial": "Liver Met., Epi_CEACAM5/6",
    # S2 E prints its two-sided P as a value through p_label, as every box
    # panel does (the page printed "p = 0.19"; the predecessor set "ns").
    "ns": "P = 0.19",
    # S2 H: each box is titled with the cohort, as the page set it; the
    # predecessor put the cohort in a figure-level suptitle (removed below).
    "CEACAM5 (Epi. deconvolved)": "TCGA-STAD\nCEACAM5 (Epi. deconvolved)",
    "CEACAM6 (Epi. deconvolved)": "TCGA-STAD\nCEACAM6 (Epi. deconvolved)",
    # S4 A-E: the peritoneal-metastasis section is labelled as the printed
    # page labels it. The predecessors set an asterisk pointing at a footnote;
    # the submitted S4 was hand-edited to drop both (00_GROUND_TRUTH/README.md,
    # amendment 2026-09-10) and the page is the ground truth.
    "GC6-PM*": "GC6-PM",
    # S2 B and C: no mathtext. A subscript or superscript glyph renders at
    # 0.7 of the size and prints under the 6 pt floor (4.2 pt on the page),
    # so log2 is set plain and a P below 0.001 is printed as the paper
    # prints it. Read old -> new.
    "Epi. CEACAM5 (log$_2$+1)": "Epi. CEACAM5 (log2+1)",
    "Epi. CEACAM6 (log$_2$+1)": "Epi. CEACAM6 (log2+1)",
    "CEACAM5 (log$_2$+1)": "CEACAM5 (log2+1)",
    "CEACAM6 (log$_2$+1)": "CEACAM6 (log2+1)",
    "ρ = 0.89, P = 2$\\times 10^{-16}$\nn = 45": "ρ = 0.89, P < 0.001\nn = 45",
    "ρ = 0.86, P = 4$\\times 10^{-121}$\nn = 410": "ρ = 0.86, P < 0.001\nn = 410",
    # S4 F, since 2026-09-16: the panel went from the page width to the right
    # half of a 3 x 2 page, and its 75 mm title breaks into two lines. The
    # words do not change.
    "GraphST Deconvolution: Cell Type Proportions per Sample":
        "GraphST Deconvolution:\nCell Type Proportions per Sample",
}

REMOVALS_S1_S6 = {
    # S1 B and C stack their three QC metrics, each row 12 mm tall; the
    # predecessor's 7 pt y labels are 18-20 mm of rotated type and cannot
    # stand beside a 12 mm row. Each row's title names the metric ("Genes
    # per Cell", "Total Counts", "Mitochondrial %"), so the y label goes.
    "Number of genes": "the row title names the metric",
    "Total UMI counts": "the row title names the metric",
    "% mitochondrial": "the row title names the metric",
    # The stacked rows share one sample axis; "Sample" stands once, under the
    # bottom row, where the predecessor printed it under each of three boxes.
    "Sample": "one shared sample axis; the label stands under the bottom row",
    # S2 F: the printed page sets no title over the metaprogram heat map.
    "Metaprogram Gene Signatures": "the printed page sets no title over S2 F",
    # S2 H: the cohort moved into each box's title (rename above).
    "TCGA-STAD: Overall Survival by Epithelial Expression (Median Split)":
        "the cohort is in each box's title",
    # S4 A-E: the footnote the printed page dropped (see "GC6-PM*" above).
    "*GC6-PM = Peritoneal metastasis (paired with GC6 primary)":
        "the printed page carries no footnote; the label is GC6-PM",
    # S1 E, since 2026-09-16 (the author's fifth reading): the redraw of
    # 2026-09-15 kept the predecessor's key, which the published S1 E never
    # had, and its second entry read "Doublets (n=0)" beside an empty
    # histogram. The key is gone; the count stands in the figure legend.
    "Singlets (n=542,121)": "the published S1 E has no key; the count is in the legend",
    "Doublets (n=0)": "the published S1 E has no key; the object holds no predicted doublets",
}

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
    ("RENAMES_S1_S6", RENAMES_S1_S6),
    ("RENAMES_SUPPLEMENTARY", RENAMES_SUPPLEMENTARY),
])

REMOVALS = _merge("REMOVALS", [
    ("REMOVALS_FIGURE_2", REMOVALS_FIGURE_2),
    ("REMOVALS_FIGURE_3", REMOVALS_FIGURE_3),
    ("REMOVALS_FIGURE_4", REMOVALS_FIGURE_4),
    ("REMOVALS_FIGURE_5", REMOVALS_FIGURE_5),
    ("REMOVALS_S1_S6", REMOVALS_S1_S6),
    ("REMOVALS_SUPPLEMENTARY", REMOVALS_SUPPLEMENTARY),
])


#: The same sections, kept apart, so a gate can ask for one figure's
#: declarations. The supplementary section is merged in wherever it is used:
#: supplementary panels are gated by the same tool as the main figures.
BY_FIGURE = {
    "1": (RENAMES_FIGURE_1, REMOVALS_FIGURE_1),
    "2": (RENAMES_FIGURE_2, REMOVALS_FIGURE_2),
    "3": (RENAMES_FIGURE_3, REMOVALS_FIGURE_3),
    "4": (RENAMES_FIGURE_4, REMOVALS_FIGURE_4),
    "5": (RENAMES_FIGURE_5, REMOVALS_FIGURE_5),
    # The supplementary gate reads the shared supplementary section and the
    # S1-S6 redraw's own; the main figures read only the shared one.
    "S": ({**RENAMES_SUPPLEMENTARY, **RENAMES_S1_S6},
          {**REMOVALS_SUPPLEMENTARY, **REMOVALS_S1_S6}),
}
