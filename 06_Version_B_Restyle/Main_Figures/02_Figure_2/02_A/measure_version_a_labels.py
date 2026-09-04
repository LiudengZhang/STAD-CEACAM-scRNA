#!/usr/bin/env python3
"""
Measure where Version A puts Figure 2's printed panel A labels, in DATA
coordinates, and write the result to `label_pins_version_a.json`.

WHY THIS EXISTS
    Version A places the nine on-plot centroid labels with `adjust_text`, which
    repels labels by their *rendered* bounding box. Version A renders on a
    280 x 280 mm canvas; Version B renders on the 66 x 68 mm box the panel
    prints in. Same algorithm, different canvas, different answer - 55 of the
    plotted positions moved and `Stem_SPINK4` landed over the wrong arm of the
    UMAP. That is a change in what the panel *shows*, so Version B cannot run
    the algorithm; it has to be told the answer.

    This script obtains that answer by measurement rather than by judgement: it
    runs Version A's own script and reads the nine label positions and the six
    leader lines back out of the figure Version A built.

HOW VERSION A IS RUN WITHOUT WRITING TO VERSION A
    Under `compare_panel_content.capture(writes="block")`, which is the same
    harness the content gate uses:
      * `Figure.savefig` is intercepted, so Version A's .png/.svg on disk are
        not rewritten;
      * every other write (`open` in a writing mode, `to_csv`, `write_text`,
        `np.save`, ...) raises.
    Version A's script is imported and executed as-is. Not one byte of
    `03_Revised_Panels/` is opened for writing. `_audit/verify_frozen.py`
    is the proof of that, run before and after.

WHAT IS MEASURED
    texts[i]   the nine labels: string + final position in data coordinates,
               after adjust_text has moved them.
    patches[i] the six leader lines adjustText drew, as the vertices of their
               path in data coordinates. Six, not nine, because adjustText only
               draws a leader for a label it actually moved far enough.
    centroids  the per-cluster median UMAP position - the position each label
               starts at, before adjust_text.

    The leader-line vertices are recorded because they are ink on the published
    page and the content gate compares them. Their endpoints are where the
    line was clipped against the label's rendered box, so they too are a
    function of Version A's canvas and cannot be recomputed at print size.

    conda run -n Liudeng_Python_310 python measure_version_a_labels.py
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
VERSION_A = (ROOT / "03_Revised_Panels" / "Main_Figures" / "02_Figure_2" /
             "02_A" / "create_epithelial_umap.py")

sys.path.insert(0, str(ROOT / "10_Reproduction"))
import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from compare_panel_content import capture                      # noqa: E402


def run_version_a():
    """Return the figure Version A builds, having written nothing."""
    import runpy
    plt.close("all")
    argv, path0 = sys.argv[:], sys.path[:]
    sys.argv = [str(VERSION_A)]
    with capture(writes="block") as cap:
        try:
            runpy.run_path(str(VERSION_A), run_name="__main__")
        finally:
            sys.argv = argv
            sys.path = path0
    figs = cap.figures or [plt.figure(n) for n in plt.get_fignums()]
    if len(figs) != 1:
        raise SystemExit(f"expected 1 figure from Version A, got {len(figs)}")
    axes = [a for a in figs[0].axes]
    if len(axes) != 1:
        raise SystemExit(f"expected 1 axes from Version A, got {len(axes)}")
    return figs[0], axes[0]


def main():
    if not VERSION_A.exists():
        raise SystemExit(f"Version A script not found: {VERSION_A}")
    print(f"running Version A under capture(writes='block'):\n  {VERSION_A}")
    fig, ax = run_version_a()

    labels = []
    for t in ax.texts:
        if isinstance(t, matplotlib.text.Annotation):
            raise SystemExit("unexpected Annotation in Version A's axes")
        x, y = t.get_position()
        labels.append({"text": t.get_text(), "x": float(x), "y": float(y)})

    leaders = []
    for p in ax.patches:
        v = np.asarray(p.get_path().vertices, dtype=float)
        leaders.append({"type": type(p).__name__,
                        "verts": [[float(a), float(b)] for a, b in v],
                        "codes": [int(c) for c in p.get_path().codes]})

    # The scatter offsets carry every cell, so the per-cluster centroid can be
    # recomputed here without loading the h5ad a second time: scanpy draws one
    # PathCollection per categorical value is NOT true (it draws one), so the
    # centroids are instead taken from the offsets grouped by facecolour.
    coll = [c for c in ax.collections if c.get_offsets() is not None
            and len(c.get_offsets())]
    offsets = np.asarray(coll[0].get_offsets(), dtype=float)
    fc = np.asarray(coll[0].get_facecolors(), dtype=float)
    print(f"  scatter: {offsets.shape[0]} points, {fc.shape[0]} facecolours")

    centroids = {}
    if fc.shape[0] == offsets.shape[0]:
        keys, inv = np.unique(np.round(fc, 6), axis=0, return_inverse=True)
        for k in range(keys.shape[0]):
            m = inv == k
            centroids["#%02x%02x%02x" % tuple(int(round(255 * c))
                                              for c in keys[k][:3])] = {
                "n": int(m.sum()),
                "median": [float(np.median(offsets[m, 0])),
                           float(np.median(offsets[m, 1]))]}

    out = {
        "source": str(VERSION_A.relative_to(ROOT)),
        "xlim": [float(v) for v in ax.get_xlim()],
        "ylim": [float(v) for v in ax.get_ylim()],
        "n_points": int(offsets.shape[0]),
        "labels": labels,
        "leaders": leaders,
        "centroids_by_colour": centroids,
        "offsets_md5_check": {
            "sum_x": float(offsets[:, 0].sum()),
            "sum_y": float(offsets[:, 1].sum())},
    }
    dest = HERE / "label_pins_version_a.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"\n{len(labels)} labels, {len(leaders)} leader lines")
    for L in labels:
        print(f"  {L['text']:<12} {L['x']:+12.7f} {L['y']:+12.7f}")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
