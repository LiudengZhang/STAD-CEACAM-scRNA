"""
Build PROVENANCE.csv: one row per printed panel of the paper.

The folder has no other way to answer "which directory holds Figure 5, panel H?".
For the main figures the directory names and the printed letters diverged after
submission and were never reconciled, and the assemblers are not a reliable
translation either: assemble_figure_2.py was re-lettered after submission and
emits A-Q where the paper prints A-N. So the printed letters come from the
figure legends, recorded here by hand in PRINTED, and the assemblers supply only
the directory and filename behind each one.

Getting this wrong has cost two retractions, both on Figure 5A. See
00_GROUND_TRUTH/README.md.

Columns
  figure                1-6, S1-S11
  printed_panel         the letter as the paper prints it
  build_path            patched | built | carried_over | schematic
  source_dir            relative to this file's directory
  source_file           the artefact the assembler consumes
  source_script         the create_*.py that writes it
  source_data           the input path constant, resolved where possible
  printed_rect_mm       x0,y0,x1,y1 in the ground-truth PDF; blank if unmeasured
  reproduces_published  yes | no | unknown | na
  note

Run: python build_provenance.py
"""

from pathlib import Path
import csv
import re
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
GROUND = ROOT / "00_GROUND_TRUTH" / "figures"
MAIN = HERE / "Main_Figures"
SUPP = HERE / "Supplementary_New"
FIXES = HERE / "Supplementary_Fixes"
# Figures 1, 4 and 6 have no panel tree in the working directory; their
# assemblers live only in the generated release.
RELEASE = ROOT / "05_Code_Release" / "github_repo" / "03_Final_Panels"

OUT = HERE / "PROVENANCE.csv"

# ---------------------------------------------------------------- printed map
# (figure, printed panel) -> assembler key(s). The assembler key equals the
# printed letter for Figures 1, 3, 4 and 5; Figure 2's assembler was re-lettered
# after submission, so its mapping is spelled out.
PRINTED = {
    "1": {"A": ["A"], "B": ["B"], "C": ["C"]},
    "2": {
        "A": ["A"],    # epithelial UMAP, 9 states
        "B": ["C"],    # CEACAM5 UMAP
        "C": ["D"],    # CEACAM6 UMAP
        "D": ["G"],    # CEACAM5 vs CEACAM6 correlation
        "E": ["H"],    # CEACAM5/6+ proportion per sample
        "F": ["B"],    # Milo
        "G": ["J"],    # metaprogram gene loading dotplot
        "H": ["K"],    # MP4 scores          (H and I share one panel file)
        "I": ["K"],    # MP5 scores
        "J": ["M"],    # checkpoint dotplot
        "K": ["N1", "N2"],
        "L": ["O1", "O2"],
        "M": ["P"],    # representative IHC
        "N": ["Q"],    # combined IHC boxplot
    },
    "3": {c: [c] for c in "ABCDEFGHIJKLMN"},
    "4": {c: [c] for c in "ABCDEFGH"},
    "5": {
        "A": ["A"], "B": ["B"], "C": ["C"], "D": ["D"], "E": ["E"],
        "F": ["F"], "G": ["G"], "H": ["H"], "I": ["I"],
        "J": ["J", "K", "L", "M"],          # CD274 in four cell types
        "K": ["N"], "L": ["O"], "M": ["P"],
        "N": ["Q1", "Q2", "Q3", "Q4"],      # GSEA dotplots, four cell types
    },
}

# Panels the assemblers do not map. assemble_figure_4.py has no entry for
# panel H, but 04_H holds it: the external validation of Mac_IL1B proportion in
# GSE239676 and GSE183904.
MANUAL = {
    ("4", "H"): ("04_H", "external_validation_boxplots.svg"),
}


def panel_swaps():
    """
    {(figure, printed letter): (directory, filename)} from the patcher.

    A panel spliced into the printed figure is described by the file that was
    spliced, not by the archived original the pre-submission assembler names.
    """
    txt = (HERE / "patch_figure_annotations.py").read_text()
    block = re.search(r"PANEL_SWAPS\s*=\s*\{(.*?)^\}", txt, re.S | re.M)
    if not block:
        return {}
    out = {}
    for fm in re.finditer(r'"Figure (\d)":\s*\[(.*?)\n    \]', block.group(1), re.S):
        fig, body = fm.group(1), fm.group(2)
        for em in re.finditer(
                r'\("([A-Z])",.*?/\s*"([^"]+)"\s*\n\s*/\s*"([^"]+)"\)', body, re.S):
            out[(fig, em.group(1))] = (em.group(2), em.group(3))
    return out


SWAPS = panel_swaps()

ASSEMBLERS = {
    "1": RELEASE / "01_Figure_1" / "assemble_figure_1.py",
    "2": MAIN / "02_Figure_2" / "assemble_figure_2.py",
    "3": MAIN / "03_Figure_3" / "assemble_figure_3.py",
    "4": RELEASE / "04_Figure_4" / "assemble_figure_4.py",
    "5": MAIN / "05_Figure_5" / "assemble_figure_5.py",
}
PANEL_ROOT = {
    "1": RELEASE / "01_Figure_1",
    "2": MAIN / "02_Figure_2",
    "3": MAIN / "03_Figure_3",
    "4": RELEASE / "04_Figure_4",
    "5": MAIN / "05_Figure_5",
}

# ------------------------------------------------- hand-recorded, per WP2
# Only panels that have actually been measured or adjudicated appear here.
# A blank rect means "not measured"; the verifier skips those rather than
# guessing at a bounding box.
RECTS = {
    ("2", "D"): "79.8,5.5,118.3,30.5",
    ("5", "A"): "2.0,3.5,65.2,45.8",
    ("5", "H"): "118.8,45.5,171.1,81.0",
}
REPRODUCES = {
    ("1", "A"): ("yes", "Replaced this round with the author's redrawn vector "
                        "package, spliced by patch_figure_annotations.py."),
    ("2", "D"): ("yes", "Replaced this round. The submitted panel's rho = 0.93 over "
                        "n = 49,696 cells came from the damaged .X: 49,696 is the "
                        "count of non-NaN rows out of 106,653 and their rho is "
                        "0.9264. From .raw: 0.44 per cell, 0.72 across kNN "
                        "metacells of 10, 0.93 across the 20 samples."),
    ("5", "A"): ("yes", "Reproduces, and the reason it did not until 2026-09-01 "
                        "was the script's input path rather than the figure. The "
                        "twelve per-cell-type *_mast_prerank_gsea.csv were moved "
                        "under GSEA/_archived/ after the figures were made and the "
                        "panel script was repointed at GSEA/post/"
                        "MoMac_gsea_hallmark.csv, a different run built on the "
                        "doubly normalised .X. Against the original table the top "
                        "nine by |NES| are the printed nine, in the printed order, "
                        "matching the published bars to 0.00014 NES. Two earlier "
                        "sessions read that disagreement the other way round and "
                        "called the published panel wrong; both were retracted."),
    ("5", "H"): ("yes", "Replaced this round, from 12_R1.8_DEG_Recompute by way of "
                        "07_R1.8_NFkB_Specificity/outputs/nfkb_per_celltype.csv - "
                        "the table Fig. S9E and S10C also read. The old source "
                        "differed in sign for B cells post (+1.01 vs -0.99) and "
                        "MoMac pre (-0.98 vs +1.06)."),
}

# Twenty-one supplementary panels have no script of their own: the analysis
# module writes the panel directly into its S*/S*_* directory. Recording them as
# script-less would leave a reader with no way back to the code, so each is
# attributed to the module that saves it. Verified by the save path in the
# module, not by a text match - grep alone points S11_A and S11_C at the wrong
# script, because two other modules mention them in a docstring.
MODULE_PANELS = {
    ("S7", "C"): "01_R1.3_Cohort_Pairing/scripts/design_defence.py",
    ("S8", "A"): "03_R1.4_MP_Direction_PrePost/scripts/mp_direction_and_prepost.py",
    ("S8", "B"): "04_R1.5_CEACAM5_vs_CEACAM6/scripts/ceacam5_vs_ceacam6.py",
    ("S8", "C"): "03_R1.4_MP_Direction_PrePost/scripts/mp_direction_and_prepost.py",
    ("S8", "D"): "03_R1.4_MP_Direction_PrePost/scripts/mp_external_validation.py",
    ("S8", "H"): "04_R1.5_CEACAM5_vs_CEACAM6/scripts/dropout_and_coexpression.py",
    ("S9", "A"): "05_R1.6_Spatial_Confounders/scripts/spatial_confounders.py",
    ("S9", "B"): "05_R1.6_Spatial_Confounders/scripts/spatial_confounders.py",
    ("S9", "C"): "06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py",
    ("S9", "D"): "06_R1.7_MoMac_Lineage_Markers/scripts/momac_lineage.py",
    ("S9", "E"): "07_R1.8_NFkB_Specificity/scripts/nfkb_specificity.py",
    ("S9", "F"): "07_R1.8_NFkB_Specificity/scripts/nfkb_specificity.py",
    ("S10", "A"): "08_R2.1_PreTx_Inflammatory/scripts/pretreatment_inflammatory.py",
    ("S10", "B"): "08_R2.1_PreTx_Inflammatory/scripts/pretreatment_inflammatory.py",
    ("S10", "C"): "08_R2.1_PreTx_Inflammatory/scripts/pretreatment_inflammatory.py",
    ("S10", "D"): "09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py",
    ("S10", "E"): "09_R2.2_Adaptive_Immune/scripts/adaptive_immune_resource.py",
    ("S11", "A"): "05_R1.6_Spatial_Confounders/scripts/spatial_positive_evidence.py",
    ("S11", "B"): "05_R1.6_Spatial_Confounders/scripts/tcga_immune_exclusion.py",
    ("S11", "C"): "07_R1.8_NFkB_Specificity/scripts/nfkb_regulon_activity.py",
    ("S11", "D"): "06_R1.7_MoMac_Lineage_Markers/scripts/celltypist_annotation.py",
}

DATA_HINT = re.compile(
    r"(?:^|[\s=(\[])((?:[A-Z][A-Z0-9_]{3,})|(?:['\"][^'\"\n]+\.(?:csv|h5ad|tsv)['\"]))")


def parse_assembler(path):
    """{assembler key: (directory, filename)} from a PANELS dict literal."""
    if not path.exists():
        sys.exit(f"assembler not found: {path}")
    out = {}
    for m in re.finditer(r"^\s*'([A-Z][0-9]?)':\s*\(\s*'([^']+)'\s*,\s*'([^']+)'",
                         path.read_text(), re.M):
        out[m.group(1)] = (m.group(2), m.group(3))
    if not out:
        sys.exit(f"no PANELS entries parsed from {path}")
    return out


def find_script(directory):
    """The create_*.py in a panel directory, or whatever single .py is there."""
    if not directory.is_dir():
        return ""
    cands = sorted(p.name for p in directory.glob("*.py")
                   if not p.name.startswith("_"))
    made = [c for c in cands if c.startswith(("create_", "run_", "make_"))]
    return "; ".join(made or cands)


def find_data(directory, script_names):
    """Path constants and literal data files the panel script reads."""
    hits = []
    for name in filter(None, script_names.split("; ")):
        p = directory / name
        if not p.exists():
            continue
        txt = p.read_text(errors="replace")
        for m in re.finditer(r"\b([A-Z][A-Z0-9_]{4,})\b", txt):
            tok = m.group(1)
            if tok.endswith(("_H5AD", "_CSV", "_DIR")) and tok not in hits:
                hits.append(tok)
        for m in re.finditer(r"['\"]([^'\"\n]*\.(?:csv|h5ad|tsv))['\"]", txt):
            v = m.group(1)
            if "output" not in v.lower() and v not in hits:
                hits.append(v)
    return "; ".join(hits[:4])


def main_figure_rows():
    rows = []
    for fig in sorted(PRINTED):
        panels = parse_assembler(ASSEMBLERS[fig])
        root = PANEL_ROOT[fig]
        for printed, keys in PRINTED[fig].items():
            dirs, files, scripts, datas = [], [], [], []
            for k in keys:
                if (fig, printed) in MANUAL:
                    d, f = MANUAL[(fig, printed)]
                elif k not in panels:
                    dirs.append(f"<unmapped:{k}>")
                    continue
                else:
                    d, f = panels[k]
                dirs.append(d)
                files.append(f"{d}/{f}")
                s = find_script(root / d)
                scripts += [f"{d}/{n}" for n in filter(None, s.split("; "))]
                datas.append(find_data(root / d, s))
            rep, note = REPRODUCES.get((fig, printed), ("unknown", ""))
            # A panel spliced into the printed figure by patch_figure_annotations.py
            # is described by the file that was spliced, not by the archived
            # original the pre-submission assembler still names.
            if (fig, printed) in SWAPS:
                sd, sf = SWAPS[(fig, printed)]
                dirs, files = [sd], [f"{sd}/{sf}"]
                sc = find_script(root / sd) if (root / sd).is_dir() else ""
                scripts = [f"{sd}/{n}" for n in filter(None, sc.split("; "))]
            if fig == "1" and printed == "A":
                dirs = ["_panel_1A"]
                files = ["_panel_1A/figure1A.pdf"]
                scripts = ["_panel_1A/build_figure.py"]
                datas = ["_panel_1A/icons_src/"]
            rows.append(dict(
                figure=fig, printed_panel=printed, build_path="patched",
                source_dir="; ".join(dict.fromkeys(dirs)),
                source_file="; ".join(dict.fromkeys(files)),
                source_script="; ".join(dict.fromkeys(filter(None, scripts))),
                source_data="; ".join(dict.fromkeys(filter(None, datas))),
                printed_rect_mm=RECTS.get((fig, printed), ""),
                reproduces_published=rep, note=note))
        # directories with no printed panel behind them
        used = {d for keys in PRINTED[fig].values() for k in keys
                if k in panels for d in [panels[k][0]]}
        used |= {d for (mf, _), (d, _f) in MANUAL.items() if mf == fig}
        for d in sorted(p.name for p in root.iterdir()
                        if p.is_dir() and not p.name.startswith("_")):
            if d in used:
                continue
            rows.append(dict(
                figure=fig, printed_panel="", build_path="patched",
                source_dir=d, source_file="",
                source_script="; ".join(
                    f"{d}/{n}"
                    for n in filter(None, find_script(root / d).split("; "))),
                source_data="", printed_rect_mm="", reproduces_published="unknown",
                note="orphan: no printed panel maps to this directory"))
    rows.append(dict(figure="6", printed_panel="", build_path="schematic",
                     source_dir="", source_file="", source_script="",
                     source_data="", printed_rect_mm="",
                     reproduces_published="na",
                     note="schematic, drawn by hand; no panel script"))
    return rows


def supplementary_rows():
    rows = []
    for s in sorted(SUPP.glob("S*/")):
        fig = s.name.split("_")[0]
        for d in sorted(p for p in s.iterdir() if p.is_dir()
                        and not p.name.startswith("_")):
            letter = d.name.split("_")[-1]
            script = find_script(d)
            if script:
                scripts = "; ".join(f"Supplementary_New/{s.name}/{d.name}/{n}"
                                    for n in filter(None, script.split("; ")))
                data = find_data(d, script)
                note = "rebuilt from this script by assemble_new_supplementaries.py"
            else:
                mod = MODULE_PANELS.get((fig, letter), "")
                scripts = f"02_New_Analyses/{mod}" if mod else ""
                data = ""
                note = ("written directly by the analysis module named in "
                        "source_script, then assembled by "
                        "assemble_new_supplementaries.py") if mod else \
                       "NO SCRIPT FOUND - attribute this panel before shipping"
            rows.append(dict(
                figure=fig, printed_panel=letter, build_path="built",
                source_dir=f"Supplementary_New/{s.name}/{d.name}",
                source_file="", source_script=scripts, source_data=data,
                printed_rect_mm="",
                reproduces_published="yes" if (script or mod) else "unknown",
                note=note))
    for name, sup in (("S1", "S1_H"),):
        d = FIXES / sup
        if d.is_dir():
            script = find_script(d)
            rows.append(dict(
                figure=name, printed_panel=sup.split("_")[-1],
                build_path="built", source_dir=f"Supplementary_Fixes/{sup}",
                source_file="",
                source_script="; ".join(
                    f"Supplementary_Fixes/{sup}/{n}"
                    for n in filter(None, script.split("; "))),
                source_data=find_data(d, script), printed_rect_mm="",
                reproduces_published="yes",
                note="the one carried-over supplementary panel that is regenerated"))
    for fig in ("S1", "S2", "S3", "S4", "S5", "S6"):
        rows.append(dict(figure=fig, printed_panel="", build_path="carried_over",
                         source_dir="00_GROUND_TRUTH/figures", source_file="",
                         source_script="", source_data="", printed_rect_mm="",
                         reproduces_published="na",
                         note="shipped unchanged from the submitted PDF"))
    return rows


def main():
    if not GROUND.is_dir():
        sys.exit(f"{GROUND} not found - create 00_GROUND_TRUTH first")
    rows = main_figure_rows() + supplementary_rows()
    derived = ["figure", "printed_panel", "build_path", "source_dir",
               "source_file", "source_script", "source_data", "printed_rect_mm",
               "note"]
    # Adjudication. Nothing in this file can regenerate any of it: the verdicts
    # took two rounds over 93 rows, some of them read by eye off the printed
    # page. reproduces_published belongs here too - the REPRODUCES map below
    # covers a handful of panels and everything else came out as "unknown".
    judged = ["reproduces_published", "verdict", "verdict_date", "evidence",
              "verdict_note"]
    cols = derived[:8] + ["reproduces_published", "note"] + judged[1:]

    def key(r):
        return (r["figure"], r["printed_panel"], r.get("source_dir", ""))

    prior = {}
    if OUT.exists():
        with open(OUT) as fh:
            for row in csv.DictReader(fh):
                prior[key(row)] = row

    def adjudicated(row):
        return any(row.get(c, "") not in ("", "unknown") for c in judged)

    for r in rows:
        was = prior.get(key(r))
        for c in judged:
            if was and was.get(c, ""):
                r[c] = was[c]            # carry the human's answer
            else:
                r.setdefault(c, "")      # a new row starts unanswered

    orphaned = [k for k, row in prior.items()
                if adjudicated(row) and k not in {key(r) for r in rows}]
    if orphaned:
        print(f"REFUSING TO WRITE {OUT.name}: {len(orphaned)} adjudicated row(s) "
              f"in the existing file have no counterpart in what this script "
              f"derives from the tree, and writing would erase them:")
        for figure, panel, sd in sorted(orphaned):
            print(f"  Figure {figure} panel {panel or '-'}  ({sd or 'no source_dir'})")
        print("Reconcile the tree with the file, or edit the file directly. "
              "This script bootstraps PROVENANCE.csv; it does not maintain it.")
        sys.exit(1)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    carried = sum(1 for r in rows if adjudicated(r))
    print(f"  carried {carried} adjudicated row(s) forward unchanged")
    print(f"{len(rows)} rows -> {OUT}")
    for tag in ("yes", "no", "unknown", "na"):
        n = sum(1 for r in rows if r["reproduces_published"] == tag)
        print(f"  reproduces_published={tag:<8} {n}")
    mism = [r for r in rows if r["build_path"] == "patched" and r["printed_panel"]
            and r["source_dir"] and not r["source_dir"].startswith("_")
            and r["source_dir"].split("; ")[0].split("_")[-1] != r["printed_panel"]]
    print(f"  printed letter != directory letter: {len(mism)} panels")
    for r in mism:
        print(f"     Figure {r['figure']} panel {r['printed_panel']:<2} "
              f"lives in {r['source_dir']}")


if __name__ == "__main__":
    main()
