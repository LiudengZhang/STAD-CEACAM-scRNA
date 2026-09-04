# The printed Figure 5A is correct. The script in this directory is not.

Read this before drawing any conclusion from `momac_pathway_enrichment.png/.svg/.pdf`
or from `create_panel_a_momac_enrichment.py`.

## What the paper prints

Figure 5A, in `00_GROUND_TRUTH/figures/Figure 5.pdf` and in the preprint
(`01_Reviewer_Materials/biorxiv_v1_708917.pdf`, page 37), shows nine Hallmark gene sets
on an axis labelled **NES (Post-NR/Post-R)**:

| enriched in non-responders (red, right) | enriched in responders (blue, left) |
|---|---|
| TNF-alpha Signaling via NF-kB (≈ +2.1)  | Coagulation |
| Inflammatory Response                    | Spermatogenesis |
| Interferon Gamma Response                | Mitotic Spindle |
| IL-6/JAK/STAT3 Signaling                 | G2-M Checkpoint |
|                                          | E2F Targets (≈ −2.6) |

This agrees with the Results sentence it supports. There is nothing wrong with it.

## What the script in this directory produces

`create_panel_a_momac_enrichment.py` reads
`Round_5/02_Preparation_for_Panels/GSEA/post/MoMac_gsea_hallmark.csv`, which the data
audit traced to the differential-expression run built on a doubly normalised `.X`
(`00_Data_Audit/FINDINGS.md`, sections 1 and 7; `MoMac.h5ad` `.X` is 45.05 % NaN).

That table gives **TNF-alpha Signaling via NF-kB as −1.77 and E2F Targets as +1.92** —
the opposite signs to the printed bars — and a different set of nine. The PNG beside it
carries yet a third set (Myc Targets, UV Response, Hypoxia) and has lost the axis
direction label. **Neither has ever been in the paper.**

## Why this file exists

Two agents in two separate sessions have looked at this directory, seen the
disagreement, and concluded that the published figure was wrong. Both conclusions were
retracted. The second attempt got as far as rebuilding the panel from a clean t-test
table before it was rolled back.

## Can it be fixed?

Not exactly, and it should not be faked. The printed nine and their signs all appear in
`../05_GSEA_Summary/gsea_data/gsea_momac.csv` — the clean Welch t-test prerank that
panel N reads, run on `.raw`. But no selection rule reproduces the printed nine: top-9
by |NES| and top-9 by nominal *P* both return Oxidative Phosphorylation, EMT, Myc
Targets V1 and Hypoxia in place of Interferon Gamma Response, IL-6/JAK/STAT3,
Coagulation, Spermatogenesis, Mitotic Spindle and G2-M Checkpoint. The magnitudes differ
too (printed ≈ 2.1 for NF-kB, 1.875 in that file). The run that produced the published
panel is not on disk.

So this panel is recorded in `../../../PROVENANCE.csv` as
`reproduces_published = no`, with the published values in the note. Anyone reproducing
the paper from the deposited code should be told that this one panel does not
regenerate, rather than shown a different panel and left to assume it does.

## If you are about to change something here

Don't, unless you have first found the original run. Changing the script to read
`gsea_momac.csv` makes it *self-consistent* but still not the published panel, and it
would silently replace four of the nine gene sets that Figure 5B and the Results
sentence both depend on.
