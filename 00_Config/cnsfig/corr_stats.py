"""
The correlation statistics a scatter panel prints, written once beside the
panel so the legend can be generated from them.

WHY THIS EXISTS (2026-09-14, evening)
    The author moved the P values of the scatter panels (Figure 3 A, B, D, E;
    Figure 5 K, M) off the panels and into the legend, and put rho inside the
    boxes. A P value that is typed into the legend is a second copy of a
    number the panel script computes (RULES.md rule 5), so the script writes
    it here - `write(panel_dir, rows)` - and 04_Manuscript_R1/01_Main_Text/
    edits.py builds the legend sentence from `sentence(figure)`.
    04_Manuscript_R1/verify_manuscript_text.py checks the legend against the
    same files, with a mutation.

    The P string is one significant figure, as the panels printed it
    ("P = 0.04", "P = 8e-7"); the printed form and the number are both kept.
"""

import csv
import math
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
#: The panel tree: Main_Figures in the working tree, 03_Final_Panels in the
#: code release (05_Code_Release/update_release.py renames it), where this
#: module is staged beside the manuscript verifier.
_MAIN = next((c for c in (_ROOT / "03_Final_Panels",
                          _ROOT / "03_Final_Panels") if c.is_dir()), None)
if _MAIN is None:
    raise FileNotFoundError(f"no panel tree beside {_ROOT}")
FILENAME = "correlation_stats.csv"

#: Printed panel -> (figure directory, panel directory), from PROVENANCE.csv
#: (the directory letter is not the printed letter: rule 2).
PANELS = {
    ("3", "A"): "03_Figure_3/03_A",
    ("3", "B"): "03_Figure_3/03_C",
    ("3", "D"): "03_Figure_3/03_B",
    ("3", "E"): "03_Figure_3/03_D",
    ("5", "K"): "05_Figure_5/05_O",
    ("5", "M"): "05_Figure_5/05_Q",
}

__all__ = ["write", "read", "p_string", "sentence", "PANELS", "FILENAME"]


def p_string(p):
    """One significant figure, floored, as the panels printed it."""
    p = float(p)
    if p <= 0:
        raise ValueError("P must be positive")
    e = math.floor(math.log10(p))
    c = int(p / 10 ** e)
    if e >= -3:
        return f"P = {c * 10 ** e:.{-e}f}"
    return f"P = {c}e{e}"


def write(panel_dir, rows):
    """rows: [(x_label, rho, p, n)] in printed order. Writes FILENAME."""
    path = Path(panel_dir) / FILENAME
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["x", "rho", "p", "p_printed", "n"])
        for x, rho, p, n in rows:
            w.writerow([x, f"{rho:.4f}", f"{p:.6g}", p_string(p), n])
    return path


def read(figure, letter):
    path = _MAIN / PANELS[(str(figure), letter)] / FILENAME
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; run the panel script")
    return list(csv.DictReader(path.open()))


def sentence(figure, letters):
    """'(A) CEACAM5, P = 0.04; CEACAM6, P = 0.1. (B) ...' from the files."""
    parts = []
    for letter in letters:
        rows = read(figure, letter)
        parts.append(f"({letter}) " + "; ".join(
            f"{r['x']}, {r['p_printed']}" for r in rows))
    return "Spearman P values: " + ". ".join(parts) + "."
