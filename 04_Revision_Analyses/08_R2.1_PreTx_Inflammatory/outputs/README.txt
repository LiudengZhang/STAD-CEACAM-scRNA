Outputs of 08_R2.1_PreTx_Inflammatory (reviewer point R2.1)

  pretx_state_abundance.csv    IL-1b+ inflammatory MoMac abundance, by group
  pretx_signature_scores.csv   that state's signature, scored per cell, by group
  il1b_signature_genes.csv     the genes that signature is built from
  nfkb_pre_vs_post.csv         TNFa/NF-kB Hallmark NES and FDR, pre and post,
                               per cell type - READ THE NOTE BELOW
  pretx_inflammatory_report.txt  the three sections written up, each with the
                               table it was rendered from named in its SOURCES

nfkb_pre_vs_post.csv IS NOT THE TABLE THE PAPER QUOTES.

  Two NF-kB enrichment runs exist and both are kept. They differ because they
  are computed on different matrices, not because one is a corrected copy of
  the other:

    live      this file. Computed by this module from
              gsea/  (12_R1.8_DEG_Recompute)
              which was produced on the .X of the input objects. Eight of those
              carry a double normalisation dated 2026-07-30 (see the data
              audit, FINDINGS.md sections 1 and 7).

    adopted   nfkb_per_celltype_sound13.csv, written by
              13_R1.8_Neutrophil_Rebuilt_Recompute. Same analysis, recomputed
              on sound input. This is the one Supplementary Table S10, the
              Results, the response letter and section 3 of the report beside
              this file all quote.

  So the numbers here will not match Supplementary Table S10, and are not meant
  to. Measured 2026-09-10 over the 52 paired values ST10's per-cell columns and
  this file have in common - 13 cell types x {pre, post} x {NES, FDR} - 43
  of the 52 differ, and MoMac pre-treatment differs in sign: here NES +1.055
  with FDR 0.472, in ST10 NES -1.001 with FDR 0.893.

  This file is kept, unaltered, because it is what this module computed and
  because one file name has to keep meaning one set of contents. Regenerating
  it from the adopted table would put a second set of contents behind a name
  that has already been cited, which is the defect this project has already
  been bitten by twice. If you want the paper's numbers, read Table S10 or
  nfkb_per_celltype_sound13.csv.
