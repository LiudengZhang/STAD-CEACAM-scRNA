#!/usr/bin/env python
"""Render CLAIM_COMPARISON.md from outputs/claim_comparison.csv.

The prose is here; every row block is generated from the CSV, so the document and
the machine-readable table cannot drift.
"""
import csv
import pathlib
from collections import Counter

HERE = pathlib.Path(__file__).resolve().parents[1]
ROWS = list(csv.DictReader((HERE / "outputs" / "claim_comparison.csv")
                           .open(encoding="utf-8")))
OUT = HERE / "CLAIM_COMPARISON.md"

ORDER = ["keep verbatim", "restate", "drop", "keep verbatim (flagged)"]
TITLE = {
    "keep verbatim": "A. Survive adoption verbatim - skim these",
    "restate": "B. Must be restated - the wording is ready to paste",
    "drop": "C. Must be dropped",
    "keep verbatim (flagged)": "D. Not covered by any measurement, or not "
                              "reproducible from any table - keep, and flag",
}


def block(r):
    L = []
    lid = r['ledger_id']
    tag = ("  ·  NOT IN ledger.csv - added here" if "not in ledger" in lid
           else (f"  ·  ledger {lid}" if lid else ""))
    L.append(f"#### {r['claim_id']}" + tag)
    L.append("")
    L.append(f"**Where** {r['location']}  ·  *{r['source']}*")
    L.append("")
    L.append("> " + r["verbatim"].replace("\n", " "))
    L.append("")
    L.append("| | |")
    L.append("|---|---|")
    fields = [
        ("As printed", r["printed_value"]),
        ("Sound-input value", r["sound_input_value"]),
        ("Table it came from", r["sound_input_table"]),
        ("**Axis 1 - input**", r["axis1_input"]),
        ("**Axis 2 - stability**", r["axis2_stability"]),
        ("Signal gate", r["gate"]),
        ("`verify_numbers.py`", r["verify_numbers_guard"]),
        ("Guard expects", r["guard_expected_value"]),
        ("Supporting measurement", r["supporting_measurement"]),
    ]
    for k, v in fields:
        if v:
            L.append(f"| {k} | {v.replace('|', '/')} |")
    L.append("")
    if r["proposed_wording"]:
        L.append("**Proposed replacement**")
        L.append("")
        L.append("> " + r["proposed_wording"].replace("\n", " "))
        L.append("")
    if r["note"]:
        L.append(r["note"])
        L.append("")
    return "\n".join(L)


def main():
    c = Counter(r["recommendation"] for r in ROWS)
    missed = [r["claim_id"] for r in ROWS if "not in ledger" in r["ledger_id"]]
    moved = [r["claim_id"] for r in ROWS if r["axis1_input"].startswith("moved")]
    guarded_moving = [r for r in ROWS
                      if r["verify_numbers_guard"].startswith("yes")
                      and "->" in r["guard_expected_value"]
                      and "-> unchanged" not in r["guard_expected_value"]]

    head = f"""# NF-κB / TNFα-signalling claims, sentence by sentence

**Module** `04_Revision_Analyses/17_NFkB_Claim_Ledger/`
**Date** 1 September 2026 · **Status** report only. Nothing is adopted, nothing is
applied. No manuscript file, figure, panel script, `PROVENANCE.csv`, `verify_numbers.py`
or sibling module's output was written; the only files created are this document, its
CSV and the two scripts that build them, all inside this directory.

**Premise, given by the author and not re-derived here.** The NF-κB analysis **adopts**
the sound-input recompute, and **Fig. S10E is to be changed**. So the question below is
not *whether* to adopt. It is: given adoption, which sentence survives, which must be
restated and to what, and which must go. Where a sentence survives unchanged that is
recorded as a finding, not as a blank row. (For the record: the on-disk file that raised
the S10E question, `03_Final_Panels/SUPPLEMENTARY_AUDIT.md`, still reads
"Needs the author's ruling. Do not repoint it unasked." It has not been updated with the
ruling; that is a bookkeeping gap, not a contradiction.)

---

## 0. The counts

| | |
|---|---|
| Claims examined | **{len(ROWS)}** |
| **Survive adoption verbatim** | **{c['keep verbatim']}** |
| **Need restating** | **{c['restate']}** |
| **Must be dropped** | **{c['drop']}** |
| **Not covered by any measurement on disk** (kept, and flagged) | **{c['keep verbatim (flagged)']}** |

Of the {len(ROWS)}, **{len(missed)} were not in `ledger.csv`** and are added here:
{', '.join(missed)}. Two of them - **C25**, the Discussion's "detectable only after
treatment", and **C29**, the Methods' subsampling sentence - need restating, so the
omission was not cosmetic.

**{len(guarded_moving)} of the claims carry a `verify_numbers.py` guard whose expected
value must move when the sentence does** - thirteen individual checks in all, listed in
section 3. Adopting a number without moving its guard turns the pre-ship gate from a
check into a guaranteed failure. A further six claims are guarded by checks that adoption
leaves passing.

A fifth count worth having: **{len(moved)} of the {len(ROWS)} claims carry a number that
actually moves** on the sound inputs. The rest either carry no number, or carry one from
a pipeline the recompute does not touch.

---

## 1. How to read a row

Every claim is scored on **two independent axes**, because they are different questions
and this project has run them together before.

**Axis 1 — does the number move when the input is sound?** Eight input h5ads carry a
double normalisation from 2025-07-30 that left whole cell rows NaN
(`00_Data_Audit/FINDINGS.md` §1 and §7). The recompute on sound per-cell-type inputs is
`04_Revision_Analyses/13_R1.8_Neutrophil_Rebuilt_Recompute/`. Two variants exist:

- **`sound12`** — twelve cell types; neutrophils have no sound source anywhere. This is
  the variant behind the recorded headline **"14 of `verify_numbers.py`'s 152 checks
  move"**, two of the fourteen coming from the missing neutrophil row and twelve from the
  matrix.
- **`sound13`** — the same twelve plus the neutrophil row carried over. **This is the
  variant used throughout below**, because it is the one `stability.csv` and the whole
  ledger are scored against, and because it keeps the denominator at thirteen so a
  "12 of 13" claim stays comparable.

The verdict is `unchanged`, `moved`, or `not covered`.

**Axis 2 — is the claim stable at all?**, independent of which input it was computed on:
seed, ranking metric, unit of replication, MAST specification, and the **signal gate** —
whether the contrast the number is read off distinguishes *any* Hallmark set from its own
permutation null. `outputs/signal.csv` scores all 52 contrast-rows on that gate;
`signal_runs.csv` holds the 365 runs behind it. The verdict is one of `safe`, `floor`,
`unstable`, `no signal`, or `gated`.

**Paragraph numbers.** Locations are given as 1-based paragraph indices over all `w:p`
elements of `Manuscript_R1_clean.docx`. `ledger.csv` numbers the same paragraphs
**0-based** — its "para 68" is this document's para 69, its "para 120" is para 121. Both
are given where it matters.

---

## 2. The two things that changed since the ledger was written

### 2.1 The S10E claim is no longer unmeasured — and it is the claim that moves most

`FINDINGS.md` §8 records Fig. S10D/E as "the one remaining unaudited GSEA claim in the
paper", and `stability.csv` scores M12 `untestable` on all five axes. **That is now out
of date.** `03_Final_Panels/SUPPLEMENTARY_AUDIT.md` fault 1 measured it:

- repointing `adaptive_immune_resource.py:146` from the prepared Round_5 tables to the
  sound recompute leaves **2 of 40 top-5 Hallmark slots** standing across the eight
  lineage × phase cells;
- `G2-M Checkpoint`, `Mitotic Spindle` and `E2F Targets` — the entire "proliferation
  programs in non-responders" reading — are in the top 5 for **neither** CD4 nor CD8 in
  **either** phase after the repoint;
- and it is not only an input problem: MAST-recompute against t-test-recompute on the
  **same sound matrix** also share **0–1 of 5**, so the lineage ranking is not stable to
  the choice of test either.

That is corrected in C07 and C41 below, and both are recommended for dropping.

### 2.2 MAST specification B is finished — and it is *not* for the manuscript

`04_Revision_Analyses/14_MAST_Specification/FINDINGS.md` (647 lines) completes the
specification `~ condition + cngeneson + (1 | sample_id)`, which converged for 26 of 26
contrasts and 102,118 of 102,130 gene fits. Under it:

| | spec B |
|---|---|
| TNFα/NF-κB positive after treatment | **12 of 13**, seven at FDR q < 0.05 |
| MoMac post | **+2.380, q < 0.001, rank 1 of 43** — the strongest of all thirteen |
| Epithelial post | +1.928, q < 0.001, rank 1 of 42 |
| B cells post | −1.665, rank 37 of 38 — still the only negative |
| MoMac **pre** | **−1.529** — the same sign flip the sound inputs give |
| agreement with the Welch t-test | 21 of 26 contrasts; **0 disagreements where either reaches q < 0.05** |

**MAST is dropped from the paper. The author has ruled this twice. None of the above is
a proposed manuscript sentence and none of it appears in any row below.** It is recorded
here as **supporting evidence held in reserve**, for one use only: if a reviewer asks why
the differential-expression method changed, spec B is the answer that a correctly
specified MAST reaches the *same* conclusion as the t-test, more strongly, and that the
deposited MAST fit was rank-deficient (13 columns, rank 12; 95.5% of genes change the
sign of the condition coefficient when the samples are renamed).

One live tension to hand to the author with it, from `VERIFICATION_ADDENDUM.md`
addenda 2 and 3: **`Manuscript_R1_clean.docx` now contains zero mentions of MAST, while
Figure 5A is drawn from a MAST prerank** (`Round_5/.../MoMac_mast_prerank_gsea.csv`,
NES +2.140, reproducing the printed panel to 0.0001) and C26 is the Results sentence that
quotes it. The figure is right — rule 1, settled 2026-09-01. The Methods no longer
describe how it was made. No remedy is proposed here; the measurements do not choose one.

---

## 3. The guards that must move with the sentences

`verify_numbers.py` is not a passive record. Every value below is an assertion that will
**fail** at the pre-ship gate the moment the analysis is adopted and the expected value is
not moved with it.

| line | check | expects | must become |
|---|---|---|---|
| :231 | `post NES with FDR < 0.05` | 4 | **3** |
| :232 | `cell types with positive pre NES` | 6 | **5** |
| :233 | `pre NES with FDR < 0.05` | 1 | **3** |
| :243 | `MoMac pre NES` (tol 0.01) | 1.06 | **−1.00** |
| :245 | `MoMac pre FDR q` (tol 0.01) | 0.47 | **0.89** |
| :247 | `Epithelial post NES` (tol 0.01) | 1.86 | **1.90** |
| :248 | `B cells post NES` (tol 0.01) | −0.99 | **−1.23** |
| :257 | `populations where TNFa/NF-kB is the top-ranked Hallmark set` | 5 | **3** |
| :261 | the named top-ranked set | {{MoMac, DC, Epithelial, Fibroblast, Mast}} | **{{MoMac, Epithelial, Fibroblast}}** |
| :264 | `Pericyte rank among Hallmark sets` | 31 | **30** |
| :266 | `B cells rank among Hallmark sets` | 37 | **38** |
| :323 | `MoMac pre NES (S10C)` | 1.06 | **−1.00** |
| :332 | S10C shows 6 of 13 pre and 12 of 13 post | 6 and 12 | **5 and 12** |

Unmoved by adoption and passing as they stand: `:228` (12 positive post), `:229`
(13 types), `:241`/`:242` (MoMac post 2.23, q 0), `:265`/`:267` (the two denominators,
49 and 38), `:324` (MoMac post S10C), `:328` (MoMac is strongest), `:275-279` (the 6–31%
epithelial cytokine range), `:345-354` (the three regulons), `:125` (NFKB1 regulon P).

**One orphan guard.** `verify_numbers.py:270` checks
`cell types where MAST and t-test agree in direction == 5`. No sentence in
`Manuscript_R1_clean.docx` or in the live response letter states it — MAST was removed
from the Methods this round and the concordance count went with it (0 occurrences of
"MAST" in the clean docx, measured). The check now guards nothing that is printed. Not
touched here; recorded for the author.

---

## 4. The claim table
"""

    parts = [head]
    for rec in ORDER:
        rows = [r for r in ROWS if r["recommendation"] == rec]
        parts.append(f"\n---\n\n## {TITLE[rec]}  ({len(rows)})\n")
        if rec == "keep verbatim":
            parts.append(
                "These are the rows the author can skim. Each one is stated to survive "
                "adoption **with no edit**; where its guard also survives untouched, "
                "that is said in the row.\n")
        if rec == "restate":
            parts.append(
                "Each carries a replacement sentence that the measurements actually "
                "support, with the measurement named. Nothing here is applied.\n")
        if rec == "drop":
            parts.append(
                "These two are the same claim, once in the manuscript and once in the "
                "letter. They must move together.\n")
        if rec == "keep verbatim (flagged)":
            parts.append(
                "**CLAUDE.md rule 1: the published figure is the ground truth and the "
                "code is the suspect.** None of these is proposed for correction. Each "
                "is recorded as not reproducible from the tables on disk, with the "
                "mechanism where one could be found and without one where it could "
                "not.\n")
        for r in rows:
            parts.append(block(r))

    parts.append(f"""
---

## 5. The positive control — the extractor was proved able to fail

A "found N sentences" from an unverified extractor is worth nothing; this project has had
four silent verification failures, one of which stood for a month. So before the sentence
set above was reported as complete, the extractor
(`work/extract_claims.py`, driven by `work/control_sentences.py`) was run three times over
the same two documents:

- **BASE** — the live text unmodified;
- **MINUS** — the same text with one known NF-κB claim sentence **deleted**. The sentence
  chosen was C01, *"Systematic GSEA across the 13 major cell types showed positive
  enrichment of the Hallmark TNFα signaling via NF-κB gene set in non-responders in 12 of
  13 populations after treatment, reaching FDR q < 0.05 in four (…; Fig. S9E)."* — the
  count the paper leads on, and the one a silent failure would most damage;
- **PLUS** — the same text with one synthetic NF-κB sentence **added**:
  *"Control sentence, not in the manuscript: the Hallmark TNF-alpha signaling via NF-kB
  set reached NES = 9.99 in an invented population."*

| source | base | on deletion | on insertion | verdict |
|---|---|---|---|---|
| `Manuscript_R1_clean.docx` | 69 sentences | **loses exactly 1**, and it is C01 | **gains exactly 1**, and it is the control sentence | **PASS** |
| response letter v3_clean | 15 sentences | n/a (C01 is not in the letter) | **gains exactly 1** | **PASS** |

Recorded in `outputs/extractor_control.csv`.

**A stated blind spot, found by the control's own logic.** The extractor keys on NF-κB /
TNFα / Hallmark / NES vocabulary. **C07 and C41 contain none of it** — "interferon-γ
programs enriched in responders and proliferation programs in non-responders" names no
NF-κB term — so the extractor does **not** find them. They are in this document because
`ledger.csv` carries them (M12, R07) and `SUPPLEMENTARY_AUDIT.md` measured them. The
sentence set above is therefore the **union** of three sources, not the output of one
tool: the extractor's 84 sentences, `ledger.csv`'s 36 rows, and the panel-level rows of
`PROVENANCE.csv` for Fig. 5A, 5H, 5I, 5N, S9E, S10C, S10E and S11C. Any future extractor
run should be judged against that union, not against itself.

---

## 6. What no restatement can fix

Unchanged from `FINDINGS.md` §10, and it bounds every row above.

1. **6 to 11 samples per contrast.** Pre-treatment contrasts have 3–4 samples per group,
   post-treatment 4–6; nineteen distinct samples in total. The per-cell Welch t-test
   treats 683 to 37,630 cells as independent replicates. Its q values describe *these two
   sets of cells*, not these two groups of patients.
2. **The pre and post sample sets are disjoint.** One patient in the whole cohort is
   sampled at both timepoints, and not into both of these sets. Every pre-versus-post
   statement in this paragraph — C13's "largely acquired on treatment" included — is a
   comparison **between different patients**. The Abstract says so; keep that sentence
   adjacent to the NF-κB paragraph, not only in the Abstract.
3. **Unequal tests on one radar.** Across Fig. 5H's 26 spokes the ranked list varies
   6.6-fold (1,024 to 6,801 genes) and the tested fraction of the 200-gene TNFα set varies
   2.2-fold (73 to 161 genes). Each NES is valid against its own null; the comparison
   *between* spokes is not like-for-like. The author has ruled that the panel is not to be
   annotated for this — so it is a fact about how the panel is **described**, and section
   B's replacement wordings are where it lands.
4. **Gene-level support is thin.** At the sample level limma-voom finds 9 genes at
   FDR < 0.05 across all 26 contrasts and DESeq2 finds 316, of which 179 are MoMac post.
   "The pathway is coordinately shifted" is supported; "these genes are differentially
   expressed" is not, outside MoMac post.
5. **Enrichment is transcriptomic, not biochemical.** Already stated, in C16 and C43.
   Those two sentences are what makes the rest defensible; they must survive every edit
   above them.

---

## 7. Files

    CLAIM_COMPARISON.md                this document
    outputs/claim_comparison.csv       the same {len(ROWS)} rows, machine-readable
    outputs/extractor_control.csv      section 5, the positive control
    scripts/build_claim_comparison.py  the row list; writes the CSV
    scripts/render_claim_comparison.py writes this document from the CSV
    work/extract_claims.py             the sentence extractor
    work/control_sentences.py          the positive control harness
    work/dump_docx.py                  paragraph dumper (read-only)

Read-only throughout: `04_Manuscript_R1/` (both .docx files, `edits.py`,
`verify_numbers.py`, the response letter), `03_Final_Panels/` (`PROVENANCE.csv`,
`SUPPLEMENTARY_AUDIT.md`), `04_Revision_Analyses/07_`, `08_`, `09_`, `12_`, `13_`, `14_`,
`15_`, `16_`, `00_Data_Audit/`, `00_GROUND_TRUTH/` and `Round_5/`. No GSEA, DEG or panel
pipeline was run. No patient name, medical record number or specimen identifier appears
in any file in this directory.
""")

    OUT.write_text("\n".join(parts), encoding="utf-8")
    print("WROTE", OUT, len("\n".join(parts).split(chr(10))), "lines")


if __name__ == "__main__":
    main()
