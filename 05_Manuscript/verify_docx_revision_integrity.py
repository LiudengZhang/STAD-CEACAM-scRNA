#!/usr/bin/env python3
"""Audit clean/tracked Word pairs, revision minimality, and embedded media.

This checker is intentionally independent of the builders.  It reads DOCX ZIP
members directly, applies or rejects Word revisions in memory, and compares the
resulting reader-visible paragraphs and drawings.  It never rewrites a DOCX.

Run from the repository root::

    python 04_Manuscript_R1/verify_docx_revision_integrity.py

Paths can be overridden for snapshots or intermediate builds; use ``--help``.
"""

from __future__ import annotations

import argparse
import hashlib
import posixpath
import re
import sys
import tempfile
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from zipfile import BadZipFile, ZipFile

from lxml import etree


W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PR = "http://schemas.openxmlformats.org/package/2006/relationships"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
Q = lambda ns, name: f"{{{ns}}}{name}"

REVISION_TAGS = {Q(W, "ins"), Q(W, "del"), Q(W, "moveFrom"), Q(W, "moveTo")}
MOVE_RANGE_TAGS = {Q(W, f"{kind}Range{edge}")
                   for kind in ("moveFrom", "moveTo")
                   for edge in ("Start", "End")}
COMMENT_TAGS = {Q(W, "commentRangeStart"), Q(W, "commentRangeEnd"),
                Q(W, "commentReference")}
EXPECTED_AUTHOR = "Liudeng Zhang"


@dataclass(frozen=True)
class State:
    paragraphs: tuple[str, ...]
    drawings: tuple[tuple[str, str, str], ...]


def _document_xml(path: Path) -> etree._Element:
    with ZipFile(path) as zf:
        return etree.fromstring(zf.read("word/document.xml"))


def _paragraph_marked(paragraph: etree._Element, tag: str) -> bool:
    ppr = paragraph.find(Q(W, "pPr"))
    if ppr is None:
        return False
    rpr = ppr.find(Q(W, "rPr"))
    return rpr is not None and rpr.find(Q(W, tag)) is not None


def _merge_forward(paragraphs: list[etree._Element]) -> None:
    """Apply Word's merge semantics for a revised paragraph mark."""
    for paragraph in paragraphs:
        parent = paragraph.getparent()
        if parent is None:
            continue
        nxt = paragraph.getnext()
        while nxt is not None and nxt.tag != Q(W, "p"):
            nxt = nxt.getnext()
        if nxt is None:
            continue
        at = 1 if nxt.find(Q(W, "pPr")) is not None else 0
        keep = [child for child in paragraph if child.tag != Q(W, "pPr")]
        for child in reversed(keep):
            nxt.insert(at, child)
        parent.remove(paragraph)


def transform_revisions(root: etree._Element, accept: bool) -> etree._Element:
    """Return a copy with all revisions accepted or rejected."""
    root = deepcopy(root)
    # A moved paragraph exists twice in the physical XML.  Acceptance keeps
    # the move destination; rejection keeps the source.  The range markers
    # bind all marked paragraphs into one named move but carry no visible text.
    discard_move = "moveFrom" if accept else "moveTo"
    keep_move = "moveTo" if accept else "moveFrom"
    moved_out = []
    # Block moves are delimited at body level.  Include every paragraph in the
    # discarded range, including newly inserted/deleted paragraphs whose mark
    # already carries w:ins/w:del and therefore cannot also carry a move mark.
    start_tag, end_tag = (Q(W, f"{discard_move}RangeStart"),
                          Q(W, f"{discard_move}RangeEnd"))
    for start in list(root.iter(start_tag)):
        rid = start.get(Q(W, "id"))
        node = start.getnext()
        while node is not None and not (node.tag == end_tag
                                        and node.get(Q(W, "id")) == rid):
            if node.tag == Q(W, "p"):
                moved_out.append(node)
            node = node.getnext()
    # Inline/legacy moves may have paragraph marks without body-level ranges.
    moved_out += [p for p in root.iter(Q(W, "p"))
                  if _paragraph_marked(p, discard_move) and p not in moved_out]
    for paragraph in moved_out:
        parent = paragraph.getparent()
        if parent is not None:
            parent.remove(paragraph)
    for node in list(root.iter()):
        if node.tag in MOVE_RANGE_TAGS and node.getparent() is not None:
            node.getparent().remove(node)
    # Inline move revisions use wrappers analogous to insertion/deletion.
    for node in list(root.iter(Q(W, discard_move))):
        parent = node.getparent()
        if parent is not None:
            parent.remove(node)
    for node in list(root.iter(Q(W, keep_move))):
        parent = node.getparent()
        if parent is None:
            continue
        at = parent.index(node)
        for child in reversed(list(node)):
            parent.insert(at, child)
        parent.remove(node)

    merge_tag = "del" if accept else "ins"
    marked = [p for p in root.iter(Q(W, "p")) if _paragraph_marked(p, merge_tag)]
    discard = Q(W, "del" if accept else "ins")
    keep = Q(W, "ins" if accept else "del")

    for node in list(root.iter(discard)):
        parent = node.getparent()
        if parent is not None:
            parent.remove(node)
    for node in list(root.iter(keep)):
        parent = node.getparent()
        if parent is None:
            continue
        at = parent.index(node)
        for child in reversed(list(node)):
            parent.insert(at, child)
        parent.remove(node)

    # Deleted text becomes ordinary text when revisions are rejected.
    if not accept:
        for text in root.iter(Q(W, "delText")):
            text.tag = Q(W, "t")
        for text in root.iter(Q(W, "delInstrText")):
            text.tag = Q(W, "instrText")

    for prop_name in ("pPr", "rPr"):
        for prop in root.iter(Q(W, prop_name)):
            for mark in list(prop):
                if mark.tag in REVISION_TAGS:
                    prop.remove(mark)
    for node in list(root.iter()):
        if node.tag in COMMENT_TAGS and node.getparent() is not None:
            node.getparent().remove(node)
    _merge_forward(marked)
    return root


def _rels_name(part: str) -> str:
    directory, base = posixpath.split(part)
    return posixpath.join(directory, "_rels", base + ".rels")


def _relationship_map(zf: ZipFile, part: str) -> dict[str, tuple[str, str]]:
    rel_name = _rels_name(part)
    if rel_name not in zf.namelist():
        return {}
    rels = etree.fromstring(zf.read(rel_name))
    out = {}
    for rel in rels:
        target = rel.get("Target", "")
        if rel.get("TargetMode") == "External":
            out[rel.get("Id")] = (target, rel.get("Type", ""))
            continue
        resolved = posixpath.normpath(posixpath.join(posixpath.dirname(part), target))
        out[rel.get("Id")] = (resolved, rel.get("Type", ""))
    return out


def _visible_paragraphs(root: etree._Element) -> tuple[str, ...]:
    rows = []
    for paragraph in root.iter(Q(W, "p")):
        text = "".join((node.text or "") for node in paragraph.iter()
                       if node.tag in (Q(W, "t"), Q(W, "tab"), Q(W, "br")))
        # Empty layout paragraphs are not stable under acceptance of a deleted
        # paragraph mark: Word may merge or retain them depending on the
        # adjacent table/section boundary.  Reader-visible paragraph content
        # and its order are stable and are the contract checked here.
        if text.strip():
            rows.append(text)
    return tuple(rows)


def state(path: Path, mode: str) -> State:
    with ZipFile(path) as zf:
        root = etree.fromstring(zf.read("word/document.xml"))
        if mode == "accept":
            root = transform_revisions(root, True)
        elif mode == "reject":
            root = transform_revisions(root, False)
        elif mode != "plain":
            raise ValueError(mode)
        rels = _relationship_map(zf, "word/document.xml")
        drawings = []
        for drawing in root.iter(Q(W, "drawing")):
            for blip in drawing.iter(Q(A, "blip")):
                rid = blip.get(Q(R, "embed")) or blip.get(Q(R, "link"))
                target, _ = rels.get(rid, ("<missing>", ""))
                digest = (hashlib.sha256(zf.read(target)).hexdigest()
                          if target in zf.namelist() else "<missing>")
                drawings.append((rid or "<missing>", target, digest))
        return State(_visible_paragraphs(root), tuple(drawings))


def _first_difference(left: tuple, right: tuple) -> str:
    for i, (a, b) in enumerate(zip(left, right)):
        if a != b:
            return f"first difference at item {i}: {a!r} != {b!r}"
    return f"different lengths: {len(left)} != {len(right)}"


def compare_states(label: str, observed: State, expected: State) -> list[str]:
    errors = []
    if observed.paragraphs != expected.paragraphs:
        errors.append(f"{label}: paragraph state differs; "
                      + _first_difference(observed.paragraphs, expected.paragraphs))
    # Relationship ids and filenames may legitimately be repacked.  Compare
    # ordered image content, which is what a reader sees.
    obs_images = tuple(x[2] for x in observed.drawings)
    exp_images = tuple(x[2] for x in expected.drawings)
    if obs_images != exp_images:
        errors.append(f"{label}: drawing state differs; "
                      + _first_difference(obs_images, exp_images))
    return errors


def _normal(text: str) -> str:
    # Collapse layout whitespace, but preserve case: changing a reference title
    # from title case to the cited article's sentence case is a real edit.
    return " ".join(text.split())


def revision_noise(path: Path) -> list[str]:
    """Find delete/insert replacements that retrack unchanged boundary text."""
    root = _document_xml(path)
    errors = []
    for p_index, paragraph in enumerate(root.iter(Q(W, "p"))):
        for parent in paragraph.iter():
            children = list(parent)
            for left, right in zip(children, children[1:]):
                if left.tag != Q(W, "del") or right.tag != Q(W, "ins"):
                    continue
                old = "".join(left.itertext())
                new = "".join(right.itertext())
                if not old and not new:
                    continue
                if _normal(old) == _normal(new):
                    errors.append(f"{path.name}: paragraph {p_index} has identical "
                                  f"delete/insert text {old!r}")
                    continue
    return errors


def author_integrity(path: Path) -> list[str]:
    """Require Liudeng Zhang in revision and core document metadata."""
    errors = []
    with ZipFile(path) as zf:
        for name in zf.namelist():
            if name.endswith(".xml") and b"Claude" in zf.read(name):
                errors.append(f"{path.name}: Claude remains in {name}")
        root = etree.fromstring(zf.read("word/document.xml"))
        authors = {node.get(Q(W, "author")) for node in root.iter()
                   if node.tag in REVISION_TAGS and node.get(Q(W, "author"))}
        if authors and authors != {EXPECTED_AUTHOR}:
            errors.append(f"{path.name}: revision authors are {sorted(authors)!r}, "
                          f"expected only {EXPECTED_AUTHOR!r}")
        core = etree.fromstring(zf.read("docProps/core.xml"))
        creator = core.find("{http://purl.org/dc/elements/1.1/}creator")
        editor = core.find("{http://schemas.openxmlformats.org/package/2006/metadata/core-properties}lastModifiedBy")
        for label, node in (("creator", creator), ("lastModifiedBy", editor)):
            value = node.text if node is not None else None
            if value != EXPECTED_AUTHOR:
                errors.append(f"{path.name}: core {label} is {value!r}, "
                              f"expected {EXPECTED_AUTHOR!r}")
    return errors


def media_integrity(path: Path) -> list[str]:
    """Check every embedded image, unused media, and duplicate media payload."""
    errors = []
    with ZipFile(path) as zf:
        names = set(zf.namelist())
        media = {name for name in names if name.startswith("word/media/")
                 and not name.endswith("/")}
        referenced = set()
        image_rels = set()
        for part in sorted(name for name in names if name.startswith("word/")
                           and name.endswith(".xml") and "/_rels/" not in name):
            try:
                root = etree.fromstring(zf.read(part))
            except etree.XMLSyntaxError as exc:
                errors.append(f"{path.name}: malformed {part}: {exc}")
                continue
            rels = _relationship_map(zf, part)
            for rid, (target, rel_type) in rels.items():
                if rel_type.endswith("/image") and not target.startswith(("http:", "https:")):
                    image_rels.add(target)
            for node in root.iter():
                for attr in (Q(R, "embed"), Q(R, "link")):
                    rid = node.get(attr)
                    if not rid:
                        continue
                    if rid not in rels:
                        errors.append(f"{path.name}: {part} uses missing relationship {rid}")
                        continue
                    target, rel_type = rels[rid]
                    if rel_type.endswith("/image"):
                        referenced.add(target)
                        if target not in names:
                            errors.append(f"{path.name}: {part} {rid} targets missing {target}")
        for target in sorted(image_rels - referenced):
            errors.append(f"{path.name}: unused image relationship targets {target}")
        for target in sorted(media - referenced):
            errors.append(f"{path.name}: orphan media member {target}")
        hashes = defaultdict(list)
        for target in sorted(media):
            hashes[hashlib.sha256(zf.read(target)).hexdigest()].append(target)
        for members in hashes.values():
            if len(members) > 1:
                errors.append(f"{path.name}: duplicate media payload: {', '.join(members)}")
    return errors


def relationship_integrity(path: Path) -> list[str]:
    """Require every internal OPC relationship to resolve to a package part."""
    errors = []
    with ZipFile(path) as zf:
        names = set(zf.namelist())
        for rel_name in sorted(name for name in names if name.endswith(".rels")):
            try:
                rels = etree.fromstring(zf.read(rel_name))
            except etree.XMLSyntaxError as exc:
                errors.append(f"{path.name}: malformed {rel_name}: {exc}")
                continue
            if rel_name == "_rels/.rels":
                base_dir = ""
            elif "/_rels/" in rel_name:
                directory, leaf = rel_name.rsplit("/_rels/", 1)
                base_part = posixpath.join(directory, leaf[:-5])
                base_dir = posixpath.dirname(base_part)
            else:
                errors.append(f"{path.name}: invalid relationship part name {rel_name}")
                continue
            for rel in rels:
                if rel.get("TargetMode") == "External":
                    continue
                target = posixpath.normpath(posixpath.join(base_dir,
                                                           rel.get("Target", "")))
                if target not in names:
                    errors.append(f"{path.name}: {rel_name} relationship "
                                  f"{rel.get('Id')} targets missing {target}")
    return errors


FORBIDDEN = (
    ("alternative response classification",
     re.compile(r"alternative[ -](?:response[ -])?classif", re.I)),
    ("P02/P26 classification or sensitivity language",
     re.compile(r"(?:P0?2|P26).{0,120}(?:classif|reclass|label|sensitiv|respon)|"
                r"(?:classif|reclass|label|sensitiv|respon).{0,120}(?:P0?2|P26)", re.I | re.S)),
    ("obsolete adjusted-density Figure S9A citation",
     re.compile(r"β\s*=\s*\+?0\.02\s*,\s*P\s*=\s*0\.68\s*;\s*Fig\.\s*S9A", re.I)),
    ("obsolete TCGA Figure S9B/S9C citation",
     re.compile(r"β\s*=\s*[−-]80\.5.{0,100}Fig\.\s*S9(?:B|B\s*(?:,|and)\s*S9C)", re.I | re.S)),
)


def forbidden_language(path: Path) -> list[str]:
    accepted = state(path, "accept")
    text = "\n".join(accepted.paragraphs)
    errors = []
    for label, pattern in FORBIDDEN:
        match = pattern.search(text)
        if match:
            sample = " ".join(match.group(0).split())[:220]
            errors.append(f"{path.name}: {label}: {sample!r}")
    return errors


def no_revisions(path: Path) -> list[str]:
    root = _document_xml(path)
    count = sum(1 for node in root.iter() if node.tag in REVISION_TAGS)
    return [f"{path.name}: clean file carries {count} revision elements"] if count else []


def no_move_markup(path: Path) -> list[str]:
    """Reject duplicated section moves in either manuscript deliverable."""
    root = _document_xml(path)
    count = sum(1 for node in root.iter()
                if node.tag in REVISION_TAGS | MOVE_RANGE_TAGS
                and etree.QName(node).localname.startswith("move"))
    return [f"{path.name}: carries {count} move-revision elements"] if count else []


def positive_controls() -> list[str]:
    """Prove the core comparisons and detectors can fail in this run."""
    convicted = []
    base = State(("alpha",), (("rId1", "word/media/a.png", "hash-a"),))
    for name, mutant in (
        ("accepted/rejected text disagreement", State(("beta",), base.drawings)),
        ("drawing disagreement", State(base.paragraphs,
                                        (("rId1", "word/media/a.png", "hash-b"),))),
    ):
        if compare_states("mutation", mutant, base):
            convicted.append(name)
    fixture = etree.fromstring(
        f'<w:p xmlns:w="{W}"><w:del><w:r><w:delText>same</w:delText></w:r></w:del>'
        f'<w:ins><w:r><w:t>same</w:t></w:r></w:ins></w:p>')
    old_loader = globals()["_document_xml"]
    try:
        globals()["_document_xml"] = lambda _path: fixture
        if revision_noise(Path("mutation.docx")):
            convicted.append("identical delete/reinsert")
        if no_revisions(Path("mutation.docx")):
            convicted.append("revision in clean file")
    finally:
        globals()["_document_xml"] = old_loader
    accepted = transform_revisions(fixture, True)
    rejected = transform_revisions(fixture, False)
    if (_visible_paragraphs(accepted) == ("same",)
            and _visible_paragraphs(rejected) == ("same",)):
        convicted.append("accept/reject transform")

    # A deliberately orphaned image proves that the package/media traversal is
    # live rather than merely counting whatever drawings happen to be present.
    with tempfile.NamedTemporaryFile(suffix=".docx") as tmp:
        with ZipFile(tmp.name, "w") as zf:
            zf.writestr("word/document.xml",
                        f'<w:document xmlns:w="{W}"><w:body/></w:document>')
            zf.writestr("word/_rels/document.xml.rels",
                        f'<Relationships xmlns="{PR}"><Relationship Id="rId1" '
                        'Type="http://schemas.openxmlformats.org/officeDocument/'
                        '2006/relationships/image" Target="media/orphan.png"/>'
                        '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/'
                        'officeDocument/2006/relationships/comments" '
                        'Target="comments-missing.xml"/>'
                        '</Relationships>')
            zf.writestr("word/media/orphan.png", b"not-a-real-image")
        if media_integrity(Path(tmp.name)):
            convicted.append("orphan media")
        if relationship_integrity(Path(tmp.name)):
            convicted.append("broken package relationship")
    if FORBIDDEN[0][1].search("under an alternative response-classification scenario"):
        convicted.append("forbidden-language pattern")
    expected = {"accepted/rejected text disagreement", "drawing disagreement",
                "identical delete/reinsert", "revision in clean file",
                "accept/reject transform", "orphan media",
                "broken package relationship",
                "forbidden-language pattern"}
    missing = expected - set(convicted)
    return ([f"positive controls did not convict: {sorted(missing)}"] if missing else [])


def defaults(root: Path) -> dict[str, Path]:
    return {
        "manuscript_clean": root / "04_Manuscript_R1/01_Main_Text/Manuscript_R1_clean.docx",
        "manuscript_tracked": root / "04_Manuscript_R1/01_Main_Text/Manuscript_R1_tracked.docx",
        "manuscript_baseline": root / "01_Reviewer_Materials/Manuscript_submitted_frozen.docx",
        "response_clean": root / "04_Manuscript_R1/05_Response_to_Reviewers/Response_to_Reviewers_CIR260753ET_v3_clean.docx",
        "response_tracked": root / "04_Manuscript_R1/05_Response_to_Reviewers/Response_to_Reviewers_CIR260753ET_v3.docx",
        "response_baseline": root / "00_Inbox/Response_to_Reviewers_CIR260753ET_v2_commented (1).docx",
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    dflt = defaults(root)
    parser = argparse.ArgumentParser(description=__doc__)
    for name, value in dflt.items():
        parser.add_argument("--" + name.replace("_", "-"), type=Path, default=value)
    args = parser.parse_args()
    paths = {name: getattr(args, name) for name in dflt}
    errors = []
    for label, path in paths.items():
        if not path.is_file():
            errors.append(f"{label}: missing {path}")
    if errors:
        print("\n".join(f"FAIL: {error}" for error in errors))
        return 1

    try:
        manuscript_clean = state(paths["manuscript_clean"], "plain")
        manuscript_accepted = state(paths["manuscript_tracked"], "accept")
        manuscript_rejected = state(paths["manuscript_tracked"], "reject")
        manuscript_baseline = state(paths["manuscript_baseline"], "plain")
        response_clean = state(paths["response_clean"], "plain")
        response_accepted = state(paths["response_tracked"], "accept")
        response_rejected = state(paths["response_tracked"], "reject")
        # The protected PI baseline is its own accepted, comment-free reading.
        response_baseline = state(paths["response_baseline"], "accept")
    except (BadZipFile, KeyError, etree.XMLSyntaxError) as exc:
        print(f"FAIL: could not parse DOCX: {exc}")
        return 1

    errors += compare_states("manuscript accept-all versus clean",
                             manuscript_accepted, manuscript_clean)
    errors += compare_states("manuscript reject-all versus submitted baseline",
                             manuscript_rejected, manuscript_baseline)
    errors += compare_states("response accept-all versus clean",
                             response_accepted, response_clean)
    errors += compare_states("response reject-all versus flattened PI baseline",
                             response_rejected, response_baseline)

    for key in ("manuscript_clean", "response_clean"):
        errors += no_revisions(paths[key])
    for key in ("manuscript_clean", "manuscript_tracked"):
        errors += no_move_markup(paths[key])
    for key in ("manuscript_tracked", "response_tracked"):
        errors += revision_noise(paths[key])
    for key in ("manuscript_clean", "manuscript_tracked",
                "response_clean", "response_tracked"):
        errors += relationship_integrity(paths[key])
        errors += media_integrity(paths[key])
        errors += forbidden_language(paths[key])
        errors += author_integrity(paths[key])
    if len(response_clean.drawings) != 14:
        errors.append("response clean: expected 14 drawings after restoring the "
                      f"TCGA figure, found {len(response_clean.drawings)}")
    errors += positive_controls()

    print(f"Manuscript: {len(manuscript_clean.paragraphs)} paragraphs, "
          f"{len(manuscript_clean.drawings)} drawings; accept/reject states checked.")
    print(f"Response: {len(response_clean.paragraphs)} paragraphs, "
          f"{len(response_clean.drawings)} drawings; accept/reject states checked.")
    print("Positive controls: accept/reject transforms, text and drawing mismatch, "
          "clean-file revisions, revision noise, orphan media, broken package "
          "relationships, and forbidden wording were challenged in memory or "
          "a temporary DOCX.")
    if errors:
        print(f"\n{len(errors)} DOCX INTEGRITY FAILURE(S):")
        for error in errors:
            print(f"  - {error}")
        return 1
    print("DOCX REVISION INTEGRITY HOLDS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
