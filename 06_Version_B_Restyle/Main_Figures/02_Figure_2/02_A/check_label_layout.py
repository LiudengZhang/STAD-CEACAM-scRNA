#!/usr/bin/env python3
"""
Confirm, by measurement on the rendered pages, that Version B's nine panel-2A
labels sit over the same clusters as Version A's and as the published figure's,
and that none of them collides at print size.

WHY NOT JUST TRUST THE CONTENT GATE
    `compare_panel_content.py` proves the two scripts hand matplotlib the same
    numbers. It says nothing about what is on the page: a label can carry the
    right coordinates and still be printed over the wrong cluster if the canvas
    changed underneath it, which is exactly what happened when `adjust_text`
    was run at print size. So this ignores the scripts and reads the rendered
    files.

    TEST 0  THE PUBLISHED PAGE           00_GROUND_TRUTH/figures/Figure 2.pdf
        The figure as submitted is the ground truth. Panel A is clipped out of
        it at its printed rect, PyMuPDF gives the nine label strings with their
        bounding boxes, and the cluster each label is printed over is read off
        the rendered pixels underneath it.

        The published panel spells its labels with an `Epi_` prefix where the
        panel script writes them bare, and places them differently. That
        divergence is pre-existing and is already adjudicated in
        PROVENANCE.csv ("reproduced", words 91%); it is reported here, never
        acted on, and nothing in PROVENANCE.csv is touched.

    TEST 1  THE TWO SVGs                 Version A's and Version B's
        Both carry the scatter as one rasterised <image> plus nine vector label
        groups. The image element's x/y/width/height and its
        `scale(1 -1) translate(0 -T)` transform describe its placement exactly,
        so any point in SVG user units becomes a pixel of that raster. For each
        label a disc of the raster centred on the label's own white box is
        classified, and the dominant cluster is the one the label is printed
        over. Done inside each file separately: no mapping between the files,
        no reference to the h5ad, no matplotlib object.

    TEST 2  THE SAME GEOMETRY, UP TO THE AXES TRANSFORM
        An axes transform is affine, so if both SVGs place the nine labels at
        the same data coordinates then one least-squares affine maps Version
        A's nine text anchors onto Version B's with no residual. A label moved
        in data space cannot hide in an affine fit of the other eight.

    TEST 3  COLLISIONS AND CONTAINMENT AT PRINT SIZE
        Version A drew on a 280 mm canvas, so its label boxes are small
        relative to the axes; Version B's are ~4.4x larger relative to the same
        data. Pinned positions could therefore collide where Version A's did
        not. The nine white boxes are measured straight out of Version B's SVG
        - the type's real rendered extent, not an estimate - and every
        overlapping pair and every box leaving the axes is reported. That is
        what decides whether a pin needs a nudge.

HOW A DOT'S COLOUR IS CLASSIFIED
    Not by nearest RGB. The dots are drawn at alpha 0.6 and overlap, and in the
    published PDF they are faint, so the same cluster appears at every tint
    from near-white to full strength; nearest-RGB then hands most faint pixels
    to whichever palette entry is palest, which is grey (MT1E) or tan
    (Stem_SPINK4), and it did - it put four published labels over the wrong
    cluster on the first run. Alpha only scales `255 - pixel`, it does not
    rotate it, so each pixel is classified by the *direction* of its ink
    against the nine palette directions. That is alpha-invariant.

    Black type has ink direction (1,1,1), which is also grey MT1E's, so on the
    published page - where the type is baked into the pixels - the dark mask is
    dilated by 6 px at 600 dpi and those pixels are dropped before classifying.
    In the SVGs the type is vector and the raster holds only dots, so no mask
    is needed there.

    conda run -n Liudeng_Python_310 python check_label_layout.py
    conda run -n Liudeng_Python_310 python check_label_layout.py --control
"""

import base64
import io
import re
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SVG_A = (ROOT / "03_Final_Panels/02_Figure_2/02_A/"
         "epithelial_umap_minor_states.svg")
SVG_B = HERE / "epithelial_umap_minor_states.svg"
CONTROL_SVG = HERE / "_superseded" / "epithelial_umap_minor_states.adjusttext.svg"
PUBLISHED_PDF = ROOT / "00_GROUND_TRUTH" / "figures" / "Figure 2.pdf"
PUBLISHED_RECT_MM = (8.2, 7.7, 52.7, 54.2)      # LETTER_MAP.csv / PROVENANCE

LABEL_COLOUR = {
    'PTMA': '#66C2A5', 'KRT19': '#FC8D62', 'CEACAM5/6': '#8DA0CB',
    'Chief_Like': '#E78AC3', 'MUC5AC': '#A6D854', 'Stem_TPX2': '#FFD92F',
    'Stem_SPINK4': '#E5C494', 'MT1E': '#B3B3B3', 'CD74': '#8DD3C7',
}
COLOUR_LABEL = {v: k for k, v in LABEL_COLOUR.items()}
PUBLISHED_NAME = {
    "Epi_PTMA": "PTMA", "Epi_KRT19": "KRT19", "Epi_CEACAM5/6": "CEACAM5/6",
    "Epi_Chief_Like": "Chief_Like", "Epi_MUC5AC": "MUC5AC",
    "Epi_TPX2": "Stem_TPX2", "Epi_SPINK4": "Stem_SPINK4",
    "Epi_MT1E": "MT1E", "Epi_CD74": "CD74",
}

PANEL_W_MM, PANEL_H_MM = 66.0, 68.0             # Version B's canvas, for TEST 3
MARGIN = dict(left=7.0, right=7.0, top=2.0, bottom=3.5)
PT_PER_MM = 72.0 / 25.4

SAMPLE_FRAC = 0.02      # sampling disc radius, as a fraction of raster width
MIN_INK = 25            # |255 - pixel| below this carries no usable direction
TYPE_DARK = 150         # a pixel darker than this on a rendered page is type
TYPE_DILATE = 6         # px, to catch the anti-aliased halo around glyphs


def _name(c):
    """A cluster name for a sampled colour. None means the disc held no dots at
    all - the label is off the cloud, which is a finding, not a gap."""
    return "(no dots)" if c is None else COLOUR_LABEL.get(c, str(c))


def _ink_directions():
    d = {}
    for hexc in LABEL_COLOUR.values():
        c = np.array([int(hexc[k:k + 2], 16) for k in (1, 3, 5)], float)
        v = 255.0 - c
        d[hexc] = v / np.linalg.norm(v)
    return list(d), np.array(list(d.values()))


DIR_HEX, DIR_MAT = _ink_directions()


def classify(px):
    """(dominant palette colour, its share, n pixels used) for an RGB array."""
    if len(px) == 0:
        return None, 0.0, 0
    ink = 255.0 - np.asarray(px, float)
    mag = np.linalg.norm(ink, axis=1)
    keep = mag > MIN_INK
    ink, mag = ink[keep], mag[keep]
    if len(ink) == 0:
        return None, 0.0, 0
    sim = (ink / mag[:, None]) @ DIR_MAT.T
    hit = np.array([DIR_HEX[i] for i in sim.argmax(axis=1)])
    vals, counts = np.unique(hit, return_counts=True)
    k = counts.argmax()
    return vals[k], counts[k] / counts.sum(), int(counts.sum())


GROW_FRACS = (SAMPLE_FRAC, 0.03, 0.04, 0.06, 0.08)


def classify_growing(img, cx, cy, mask=None):
    """As `classify`, but widen the disc until it holds dots.

    A label need not be printed on top of its cluster to name it: the published
    panel puts "Epi_SPINK4" in the white space just off the tan arm, so a disc
    of the standard radius holds no dots at all. Reporting that as "(no dots)"
    would say nothing about which cluster the label belongs to, which is the
    question. So the disc grows until it finds ink, and the radius it needed is
    reported, so a grown sample is never mistaken for a tight one.
    """
    for frac in GROW_FRACS:
        c, share, n = classify(sample_disc(img, cx, cy, frac, mask))
        if c is not None:
            return c, share, n, frac
    return None, 0.0, 0, GROW_FRACS[-1]


def sample_disc(img, cx, cy, frac=SAMPLE_FRAC, mask=None):
    """The pixels of `img` inside a disc, minus anything `mask` marks out."""
    H, W = img.shape[:2]
    rad = max(4.0, frac * W)
    rr = np.arange(max(0, int(cy - rad)), min(H, int(cy + rad) + 1))
    cc = np.arange(max(0, int(cx - rad)), min(W, int(cx + rad) + 1))
    if not len(rr) or not len(cc):
        return np.empty((0, 3))
    sel = ((rr[:, None] - cy) ** 2 + (cc[None, :] - cx) ** 2) <= rad ** 2
    if mask is not None:
        sel = sel & ~mask[np.ix_(rr, cc)]
    return img[np.ix_(rr, cc)][sel].astype(float)


# --------------------------------------------------------------------------
# SVG
# --------------------------------------------------------------------------
def read_svg(path):
    s = Path(path).read_text()
    i = s.find("<image")
    img_tag = s[i:s.find(">", i) + 1]

    def attr(name):
        m = re.search(rf'{name}="([^"]+)"', img_tag)
        if m is None:
            raise SystemExit(f"{path}: <image> has no {name}")
        return float(m.group(1))

    b64 = re.search(r'base64,([^"]+)"', img_tag).group(1)
    raster = np.asarray(Image.open(io.BytesIO(base64.b64decode(b64)))
                        .convert("RGB"))
    T = float(re.search(r'translate\(0 (-?[\d.]+)\)', img_tag).group(1))
    image = dict(x=attr("x"), y=attr("y"), w=attr("width"), h=attr("height"),
                 T=-T, raster=raster)

    # Split on the group markers: a matplotlib text group is
    # <g id="text_N"><g id="patch_M"><path/></g><text/></g>, and a naive
    # non-greedy </g>\s*</g> swallows all nine groups in one match.
    labels = []
    for blk in [b for b in re.split(r'(?=<g id="text_\d+">)', s)
                if b.startswith('<g id="text_')]:
        m = re.search(r'<text[^>]*x="([\d.-]+)"[^>]*y="([\d.-]+)"[^>]*>'
                      r'(.*?)</text>', blk, re.S)
        if m is None:
            continue
        d = re.search(r'<path d="([^"]+)"', blk)
        verts = [(float(a), float(b)) for a, b in
                 re.findall(r'([-\d.]+)\s+([-\d.]+)', d.group(1))] if d else []
        labels.append(dict(
            text=re.sub(r"<[^>]+>", "", m.group(3)).strip(),
            anchor=(float(m.group(1)), float(m.group(2))),
            box=(min(v[0] for v in verts), min(v[1] for v in verts),
                 max(v[0] for v in verts), max(v[1] for v in verts))
            if verts else None))
    if len(labels) != 9:
        raise SystemExit(f"{path}: found {len(labels)} labels, expected 9")
    return image, labels


def svg_cluster_under(image, box):
    """Which cluster's dots lie under a label, inside one SVG."""
    H, W = image["raster"].shape[:2]
    X, Y = (box[0] + box[2]) / 2.0, (box[1] + box[3]) / 2.0
    c = (X - image["x"]) / image["w"] * W - 0.5
    r = (image["T"] - image["y"] - Y) / image["h"] * H - 0.5
    return classify(sample_disc(image["raster"], c, r))


# --------------------------------------------------------------------------
# The published page
# --------------------------------------------------------------------------
def published_panel(dpi=600):
    """(rendered crop, type mask, {label -> bbox centre in crop pixels})."""
    import fitz
    mm = 72.0 / 25.4
    doc = fitz.open(PUBLISHED_PDF)
    page = doc[0]
    x0, y0, x1, y1 = PUBLISHED_RECT_MM
    rect = fitz.Rect(x0 * mm, y0 * mm, x1 * mm, y1 * mm)
    pix = page.get_pixmap(clip=rect, dpi=dpi)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(
        pix.height, pix.width, pix.n)[:, :, :3]
    sx, sy = pix.width / rect.width, pix.height / rect.height
    found = {}
    for blk in page.get_text("dict", clip=rect)["blocks"]:
        for line in blk.get("lines", []):
            txt = "".join(sp["text"] for sp in line["spans"]).strip()
            if txt not in PUBLISHED_NAME:
                continue
            bx = [sp["bbox"] for sp in line["spans"]]
            found[PUBLISHED_NAME[txt]] = (
                ((min(b[0] for b in bx) + max(b[2] for b in bx)) / 2.0
                 - rect.x0) * sx,
                ((min(b[1] for b in bx) + max(b[3] for b in bx)) / 2.0
                 - rect.y0) * sy)
    doc.close()
    mask = ndimage.binary_dilation(np.all(img < TYPE_DARK, axis=2),
                                   iterations=TYPE_DILATE)
    missing = set(PUBLISHED_NAME.values()) - set(found)
    if missing:
        raise SystemExit(f"published panel A: could not find {sorted(missing)}")
    return img, mask, found


# --------------------------------------------------------------------------
def affine_residual(A, B):
    """Least-squares affine A->B over the nine anchors; residual in B's units."""
    P = np.array([L["anchor"] for L in A])
    Q = np.array([L["anchor"] for L in B])
    M = np.hstack([P, np.ones((len(P), 1))])
    sol, *_ = np.linalg.lstsq(M, Q, rcond=None)
    res = np.linalg.norm(M @ sol - Q, axis=1)
    return res, float(np.sqrt(abs(np.linalg.det(sol[:2, :2])))), sol


def main():
    for p in (SVG_A, SVG_B, PUBLISHED_PDF):
        if not Path(p).exists():
            raise SystemExit(f"missing {p}")
    imA, labA = read_svg(SVG_A)
    imB, labB = read_svg(SVG_B)
    byA = {L["text"]: L for L in labA}
    byB = {L["text"]: L for L in labB}
    if set(byA) != set(byB) != set(LABEL_COLOUR):
        raise SystemExit(f"label strings differ:\n  A {sorted(byA)}\n"
                         f"  B {sorted(byB)}")
    print(f"Version A svg  raster {imA['raster'].shape[1]}x"
          f"{imA['raster'].shape[0]} px, 9 labels")
    print(f"Version B svg  raster {imB['raster'].shape[1]}x"
          f"{imB['raster'].shape[0]} px, 9 labels")
    fails = []

    # ---- TEST 0 ----------------------------------------------------------
    pimg, pmask, pfound = published_panel()
    print(f"\nTEST 0  the published page: {PUBLISHED_PDF.name} panel A, "
          f"printed rect {PUBLISHED_RECT_MM} mm, {pimg.shape[1]}x"
          f"{pimg.shape[0]} px")
    print(f"{'label':<13} {'published: under':<17} {'share':>6}  "
          f"{'B: under':<13} {'share':>6}  agree")
    # "under" here means the nearest ink the label sits on or beside - see
    # classify_growing; the published panel places one label off its cluster.
    for name in LABEL_COLOUR:
        cp, sp, _, fr = classify_growing(pimg, *pfound[name], mask=pmask)
        cb, sb, _ = svg_cluster_under(imB, byB[name]["box"])
        agree = (cp == cb)
        grown = "" if fr == SAMPLE_FRAC else f"  (disc grown to {fr:g})"
        print(f"{name:<13} {_name(cp):<17} {sp:6.2f}  {_name(cb):<13} "
              f"{sb:6.2f}  {'yes' if agree else 'NO'}{grown}")
        if not agree:
            fails.append(f"TEST 0 {name}: over {_name(cp)} on the published "
                         f"page but over {_name(cb)} in Version B")
        if cp != LABEL_COLOUR[name]:
            fails.append(f"TEST 0 {name}: the PUBLISHED page prints it over "
                         f"{_name(cp)}, not its own cluster")

    # ---- TEST 1 ----------------------------------------------------------
    print("\nTEST 1  the cluster under each label, read off each SVG separately")
    print(f"{'label':<13} {'A: under':<13} {'share':>6}  "
          f"{'B: under':<13} {'share':>6}  same  own")
    for name in LABEL_COLOUR:
        ca, sa, _ = svg_cluster_under(imA, byA[name]["box"])
        cb, sb, _ = svg_cluster_under(imB, byB[name]["box"])
        same, own = (ca == cb), (cb == LABEL_COLOUR[name])
        print(f"{name:<13} {_name(ca):<13} {sa:6.2f}  {_name(cb):<13} "
              f"{sb:6.2f}  {'yes' if same else 'NO':<4}  "
              f"{'yes' if own else 'no'}")
        if not same:
            fails.append(f"TEST 1 {name}: over {_name(ca)} in Version A but "
                         f"over {_name(cb)} in Version B")
        if not own:
            fails.append(f"TEST 1 {name}: printed over {_name(cb)}, not its "
                         f"own cluster")

    # ---- TEST 2 ----------------------------------------------------------
    res, scale, _ = affine_residual([byA[n] for n in LABEL_COLOUR],
                                    [byB[n] for n in LABEL_COLOUR])
    print(f"\nTEST 2  one affine maps A's nine anchors onto B's; "
          f"scale {scale:.6f}")
    print(f"        residual  max {res.max():.6f} pt, "
          f"rms {np.sqrt((res ** 2).mean()):.6f} pt")
    if res.max() > 0.05:
        fails.append(f"TEST 2: affine residual {res.max():.4f} pt - the two "
                     f"files do not place the nine labels at the same data "
                     f"coordinates")
        for name, r in zip(LABEL_COLOUR, res):
            print(f"          {name:<13} {r:8.4f} pt")

    # ---- TEST 3 ----------------------------------------------------------
    ax = (MARGIN["left"] * PT_PER_MM, MARGIN["top"] * PT_PER_MM,
          (PANEL_W_MM - MARGIN["right"]) * PT_PER_MM,
          (PANEL_H_MM - MARGIN["bottom"]) * PT_PER_MM)
    print(f"\nTEST 3  Version B label boxes at print size, against the axes "
          f"rect ({ax[0]:.1f}, {ax[1]:.1f}) - ({ax[2]:.1f}, {ax[3]:.1f}) pt")
    names = list(LABEL_COLOUR)
    boxes = {n: byB[n]["box"] for n in names}
    print(f"{'label':<13} {'x0':>7} {'y0':>7} {'x1':>7} {'y1':>7}   "
          f"{'w mm':>6} {'h mm':>6}  inside")
    for n in names:
        b = boxes[n]
        inside = (b[0] >= ax[0] - 1e-6 and b[1] >= ax[1] - 1e-6
                  and b[2] <= ax[2] + 1e-6 and b[3] <= ax[3] + 1e-6)
        print(f"{n:<13} {b[0]:7.2f} {b[1]:7.2f} {b[2]:7.2f} {b[3]:7.2f}   "
              f"{(b[2]-b[0])/PT_PER_MM:6.2f} {(b[3]-b[1])/PT_PER_MM:6.2f}  "
              f"{'yes' if inside else 'NO'}")
        if not inside:
            fails.append(f"TEST 3 {n}: label box leaves the axes")
    hits, gaps = 0, []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = boxes[names[i]], boxes[names[j]]
            ox = min(a[2], b[2]) - max(a[0], b[0])
            oy = min(a[3], b[3]) - max(a[1], b[1])
            if ox > 0 and oy > 0:
                hits += 1
                print(f"  OVERLAP {names[i]} x {names[j]}: {ox:.2f} x {oy:.2f} pt")
                fails.append(f"TEST 3: {names[i]} overlaps {names[j]}")
            gaps.append((max(max(a[0], b[0]) - min(a[2], b[2]),
                             max(a[1], b[1]) - min(a[3], b[3])),
                         names[i], names[j]))
    gaps.sort()
    print(f"  {hits} overlapping pair(s) of 36. Tightest three clearances:")
    for g, n1, n2 in gaps[:3]:
        print(f"    {n1} / {n2}: {g:.2f} pt = {g / PT_PER_MM:.2f} mm")

    print()
    if fails:
        print(f"{len(fails)} FAILURE(S):")
        for f in fails:
            print(f"  {f}")
        return 1
    print("All four tests pass. Every label is printed over the cluster it "
          "names, on the\npublished page and in both versions; the two SVGs "
          "place the nine labels at the same\ndata coordinates; and no label "
          "box collides or leaves the axes at print size, so no\npin needed a "
          "nudge.")
    return 0


# --------------------------------------------------------------------------
def control():
    """Prove the tests can report a failure before a pass is believed.

    The specimen is real, not synthetic: `_superseded/..adjusttext.svg` is the
    Version B panel as it was drawn on 2026-09-01 with `adjust_text` still in
    it, at this same 66 x 68 mm size. That render is the defect this work
    exists to remove. If TEST 1 cannot see it, it cannot see anything, and its
    clean verdict on the pinned panel is worthless.

    WHICH label is misplaced there was itself settled by measurement, and the
    answer is not the one the stage-3 note recorded. Mapping that render's nine
    text anchors into data coordinates - its canvas is identical to the pinned
    panel's, so the pinned panel calibrates the map to 7e-8 data units:

        CEACAM5/6    moved 4.37 units; its nearest centroid becomes
                     Stem_SPINK4's (2.90) rather than its own (3.58), and no
                     dots at all lie beneath it
        MT1E         moved 2.97 units; prints over the orange KRT19 body
        MUC5AC       moved 1.80 units; prints over the pink Chief_Like body
        KRT19        moved 1.47 units; nearest centroid becomes Chief_Like's
        Stem_SPINK4  moved 0.93 units and stays nearest its OWN cluster

    So it is CEACAM5/6 that lands over the wrong arm, not Stem_SPINK4. The
    decision to stop was right; the label it was attributed to was not.

    Required: TEST 1 reports at least one label over a different cluster from
    Version A, CEACAM5/6 among them; TEST 2's residual blows up, because no one
    affine maps A's nine anchors onto these nine; and TEST 0's published-page
    sampler changes its answer when it samples somewhere else.
    """
    if not CONTROL_SVG.exists():
        raise SystemExit(f"control specimen missing: {CONTROL_SVG}")
    imA, labA = read_svg(SVG_A)
    imC, labC = read_svg(CONTROL_SVG)
    byA = {L["text"]: L for L in labA}
    byC = {L["text"]: L for L in labC}
    print(f"POSITIVE CONTROL - the pre-pin adjustText render\n  {CONTROL_SVG}")
    moved = []
    for name in LABEL_COLOUR:
        ca, _, _ = svg_cluster_under(imA, byA[name]["box"])
        cc, sc_, _ = svg_cluster_under(imC, byC[name]["box"])
        flag = "" if ca == cc else "   <-- DETECTED, differs from Version A"
        print(f"  {name:<13} A: {_name(ca):<13} control: {_name(cc):<13} "
              f"{sc_:5.2f}{flag}")
        if ca != cc:
            moved.append(name)
    res, _, _ = affine_residual([byA[n] for n in LABEL_COLOUR],
                                [byC[n] for n in LABEL_COLOUR])
    print(f"  TEST 2 residual on the control: max {res.max():.4f} pt "
          f"(the pinned panel scores 0.000001)")

    # TEST 0's sampler must depend on WHERE it samples: read the published page
    # at each label's position but with the assignments rotated by one.
    pimg, pmask, pfound = published_panel()
    order = list(pfound)
    wrong = 0
    for i, n in enumerate(order):
        c0, _, _, _ = classify_growing(pimg, *pfound[n], mask=pmask)
        c1, _, _, _ = classify_growing(
            pimg, *pfound[order[(i + 1) % len(order)]], mask=pmask)
        wrong += (c0 != c1)
    print(f"  TEST 0 sampled at rotated positions: {wrong}/9 labels report a "
          f"different cluster")

    bad = []
    if not moved:
        bad.append("TEST 1 did NOT detect any label over a different cluster")
    elif "CEACAM5/6" not in moved:
        bad.append(f"TEST 1 detected {moved} but NOT CEACAM5/6, the label "
                   f"measured to have crossed onto another cluster here")
    if res.max() < 0.05:
        bad.append(f"TEST 2 residual {res.max():.4f} pt - the affine fit did "
                   f"not notice nine moved labels")
    if wrong < 6:
        bad.append(f"TEST 0: only {wrong}/9 rotated samples changed answer - "
                   f"the published-page sampler barely depends on position")
    print()
    if bad:
        print("POSITIVE CONTROL FAILED:")
        for b in bad:
            print(f"  {b}")
        return 1
    print(f"Positive control passed. TEST 1 detects {len(moved)} label(s) over "
          f"a different cluster\n(CEACAM5/6 among them), TEST 2's residual "
          f"blows up, and TEST 0's sampler changes its\nanswer when it samples "
          f"elsewhere. A clean verdict from these tests means something.")
    return 0


if __name__ == "__main__":
    sys.exit(control() if "--control" in sys.argv else main())
