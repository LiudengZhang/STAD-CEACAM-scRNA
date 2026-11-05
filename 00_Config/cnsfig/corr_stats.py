"""
Spearman statistics for the scatter panels (Figure 3 A, B, D, E; Figure 5 K, M).

Each panel script calls write(panel_dir, rows) to save correlation_stats.csv
beside the panel; sentence(figure, letters) builds the legend's P-value
sentence from those files.
"""

import csv
import math
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
# The main-figure panel tree.
_MAIN = next((c for c in (_ROOT / "03_Final_Panels",
                          _ROOT / "03_Final_Panels") if c.is_dir()), None)
if _MAIN is None:
    raise FileNotFoundError(f"no panel tree beside {_ROOT}")
FILENAME = "correlation_stats.csv"

# Printed panel -> panel directory (directory letters differ from printed ones).
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
    """P string with one significant figure, floored rather than rounded."""
    p = float(p)
    if p <= 0:
        raise ValueError("P must be positive")
    e = math.floor(math.log10(p))
    c = int(p / 10 ** e)
    if e >= -3:
        return f"P = {c * 10 ** e:.{-e}f}"
    # Below 0.001: Unicode superscript exponent. Panels draw it through
    # cnsfig.superscript instead, since the panel font lacks a superscript minus.
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"P = {c} × 10{str(e).translate(sup)}"


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
