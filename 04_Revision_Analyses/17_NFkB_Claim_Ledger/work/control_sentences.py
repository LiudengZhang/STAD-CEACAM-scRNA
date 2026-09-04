# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
"""Positive control for the CLAIM_COMPARISON.md sentence extractor.

The extractor is run three times:
  BASE     the live text (clean manuscript + live response letter), unmodified
  MINUS    the same text with one known NF-kB claim sentence DELETED
  PLUS     the same text with one synthetic NF-kB claim sentence ADDED

The control passes only if the extractor's sentence set LOSES exactly the deleted
sentence in MINUS and GAINS exactly the added sentence in PLUS. A tool that cannot
lose a sentence cannot be trusted when it says it found them all.
"""
import sys, csv, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from extract_claims import docx_paras, scan_paras

ROOT = pathlib.Path("/path/to/Project_4_05232025/"
                    "Round_7_major_revision")
CLEAN = ROOT / "04_Manuscript_R1/01_Main_Text/Manuscript_R1_clean.docx"
LETTER = (ROOT / "04_Manuscript_R1/05_Response_to_Reviewers/"
                 "Response_to_Reviewers_CIR260753ET_v3_clean.docx")

# The sentence deliberately removed. It is claim M01/M02, the count the paper leads
# on, chosen because it is the one a silent extractor failure would most damage.
DELETE = ("Systematic GSEA across the 13 major cell types showed positive enrichment "
          "of the Hallmark TNF")
# The sentence deliberately added. It is not in either document.
ADD = ("Control sentence, not in the manuscript: the Hallmark TNF-alpha signaling via "
       "NF-kB set reached NES = 9.99 in an invented population.")


def variant(paras, mode):
    out = []
    for p in paras:
        if mode == "minus" and DELETE in p:
            # drop only the one sentence, keep the rest of the paragraph
            i = p.index(DELETE)
            j = p.find("Fig. S9E).", i)
            p = p[:i] + p[j + len("Fig. S9E)."):] if j > 0 else p[:i]
        out.append(p)
    if mode == "plus":
        out.append(ADD)
    return out


def hits(paras, label):
    return {h[3] for h in scan_paras(paras, label)}


def main():
    rows = []
    for name, path in (("main_text_clean", CLEAN), ("response_letter_v3_clean", LETTER)):
        base_p = docx_paras(path)
        base = hits(base_p, name)
        minus = hits(variant(base_p, "minus"), name)
        plus = hits(variant(base_p, "plus"), name)
        lost = base - minus
        gained = plus - base
        rows.append(dict(
            source=name,
            n_sentences_base=len(base),
            n_sentences_minus=len(minus),
            n_sentences_plus=len(plus),
            lost_on_deletion=len(lost),
            gained_on_insertion=len(gained),
            lost_sentence=("; ".join(sorted(lost))[:180] if lost else ""),
            gained_sentence=("; ".join(sorted(gained))[:180] if gained else ""),
        ))
    # the deleted sentence only exists in the main text; the added one in both
    ok_main = (rows[0]["lost_on_deletion"] == 1 and rows[0]["gained_on_insertion"] == 1)
    ok_letter = (rows[1]["gained_on_insertion"] == 1)
    for r in rows:
        r["verdict"] = "PASS" if (r["source"] == "main_text_clean" and ok_main) or \
                                 (r["source"] != "main_text_clean" and ok_letter) else "FAIL"
    out = pathlib.Path(__file__).parent.parent / "outputs" / "extractor_control.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)
    print("WROTE", out)


if __name__ == "__main__":
    main()
