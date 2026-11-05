NF-kB Hallmark summary files
============================

nfkb_rankings_13types_ttest.csv
    Welch t-test analysis (scanpy rank_genes_groups,
    method='t-test_overestim_var'). This is the analysis the Methods
    describe for Figures 5H and 5I, and the file those panels read.

nfkb_rankings_13types_mast.csv
    MAST hurdle-model analysis, used as an independent check in the
    revision. Matches the per-cell-type tables in pre/ and post/.

In the original working tree the file named '..._mast.csv' was byte-identical to the t-test file, i.e. the name
misdescribed its contents. Figures 5H/I were unaffected because the
t-test analysis is what the Methods specify for them, but the naming
was corrected here before deposition.

Sign convention for the per-cell-type tables in pre/ and post/:
    run_gsea_from_mast.py negates the MAST log fold change before
    ranking. MAST logfoldchanges is already non-responder-relative
    (verified against IL1B, IL1A, IL6, TNF, CXCL8 and NFKBIA, all of
    which are higher in post-treatment non-responders), so in the
    saved tables a NEGATIVE NES means enriched in non-responders.
    The revision analyses negate it back; see
    04_Revision_Analyses/07_R1.8_NFkB_Specificity/scripts/.
