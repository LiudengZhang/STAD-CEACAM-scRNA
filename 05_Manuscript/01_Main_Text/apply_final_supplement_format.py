#!/usr/bin/env python3
"""Keep only supplementary-figure titles in the manuscript.

The full legends are maintained in supplementary_legends.json and printed in
the combined supplementary-figure PDF.  In the tracked manuscript, removal of
each description is attributed to ``new editor 1``.
"""

from __future__ import annotations

import json
import shutil
import tempfile
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree

HERE = Path(__file__).resolve().parent
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
CP = "http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
Q = lambda ns, name: f"{{{ns}}}{name}"
AUTHOR = "new editor 1"


def legends() -> dict[str, str]:
    return json.loads((HERE / "supplementary_legends.json").read_text())


def title_of(legend: str) -> str:
    marker = legend.find(" (A")
    if marker < 0:
        raise ValueError(f"legend has no panel marker: {legend}")
    return legend[:marker]


def visible_text(paragraph) -> str:
    out = []
    for node in paragraph.iter(Q(W, "t")):
        if any(a.tag == Q(W, "del") for a in node.iterancestors()):
            continue
        out.append(node.text or "")
    return "".join(out)


def next_revision_id(root) -> int:
    ids = []
    for node in root.iter():
        value = node.get(Q(W, "id"))
        if value and value.isdigit():
            ids.append(int(value))
    return max(ids, default=0) + 1


def deletion(run, text: str, revision_id: int):
    wrapper = etree.Element(Q(W, "del"))
    wrapper.set(Q(W, "id"), str(revision_id))
    wrapper.set(Q(W, "author"), AUTHOR)
    wrapper.set(Q(W, "date"), datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"))
    copied = deepcopy(run)
    for node in copied.iter(Q(W, "t")):
        node.tag = Q(W, "delText")
    texts = list(copied.iter(Q(W, "delText")))
    if not texts:
        raise ValueError("run has no text")
    texts[0].text = text
    for node in texts[1:]:
        node.text = ""
    wrapper.append(copied)
    return wrapper


def truncate_tracked(root, wanted: dict[str, str]) -> int:
    revision_id = next_revision_id(root)
    edited = 0
    for paragraph in root.iter(Q(W, "p")):
        current = visible_text(paragraph)
        key = next((k for k in wanted if current.startswith(f"Figure {k}.")), None)
        if key is None:
            continue
        title = title_of(wanted[key])
        if current == title:
            continue
        if current != wanted[key] and not current.startswith(title + " "):
            raise ValueError(f"unexpected {key} paragraph: {current}")
        offset = 0
        for node in list(paragraph.iter(Q(W, "t"))):
            if any(a.tag == Q(W, "del") for a in node.iterancestors()):
                continue
            text = node.text or ""
            start, end = offset, offset + len(text)
            offset = end
            if end <= len(title):
                continue
            cut = max(0, len(title) - start)
            suffix = text[cut:]
            if not suffix:
                continue
            run = node.getparent()
            parent = run.getparent()
            if cut:
                node.text = text[:cut]
                insert_at = parent.index(run) + 1
            else:
                insert_at = parent.index(run)
                parent.remove(run)
            parent.insert(insert_at, deletion(run, suffix, revision_id))
            revision_id += 1
        if visible_text(paragraph) != title:
            raise ValueError(f"failed to truncate {key}: {visible_text(paragraph)!r}")
        edited += 1
    if edited not in (0, len(wanted)):
        raise ValueError(f"edited {edited} supplementary legends; expected 0 or {len(wanted)}")
    return edited


def truncate_clean(root, wanted: dict[str, str]) -> int:
    edited = 0
    for paragraph in root.iter(Q(W, "p")):
        current = "".join(n.text or "" for n in paragraph.iter(Q(W, "t")))
        key = next((k for k in wanted if current.startswith(f"Figure {k}.")), None)
        if key is None:
            continue
        title = title_of(wanted[key])
        if current == title:
            continue
        if current != wanted[key] and not current.startswith(title + " "):
            raise ValueError(f"unexpected clean {key} paragraph: {current}")
        remaining = len(title)
        for node in list(paragraph.iter(Q(W, "t"))):
            text = node.text or ""
            if remaining >= len(text):
                remaining -= len(text)
            elif remaining > 0:
                node.text = text[:remaining]
                remaining = 0
            else:
                node.text = ""
        edited += 1
    if edited not in (0, len(wanted)):
        raise ValueError(f"edited {edited} supplementary legends; expected 0 or {len(wanted)}")
    return edited


def rewrite_docx(path: Path, tracked: bool, wanted: dict[str, str]) -> int:
    with ZipFile(path) as source, tempfile.NamedTemporaryFile(suffix=".docx", delete=False) as tmp:
        temp_path = Path(tmp.name)
        with ZipFile(temp_path, "w", ZIP_DEFLATED) as target:
            edited = 0
            for info in source.infolist():
                data = source.read(info.filename)
                if info.filename == "word/document.xml":
                    root = etree.fromstring(data)
                    edited = truncate_tracked(root, wanted) if tracked else truncate_clean(root, wanted)
                    data = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone="yes")
                elif info.filename == "docProps/core.xml":
                    root = etree.fromstring(data)
                    node = root.find(Q(CP, "lastModifiedBy"))
                    if node is None:
                        node = etree.SubElement(root, Q(CP, "lastModifiedBy"))
                    node.text = AUTHOR
                    data = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone="yes")
                target.writestr(info, data)
    shutil.move(temp_path, path)
    return edited


def main() -> None:
    wanted = legends()
    clean = HERE / "Manuscript_R1_clean.docx"
    tracked = HERE / "Manuscript_R1_tracked.docx"
    print(f"clean:   {rewrite_docx(clean, False, wanted)} descriptions removed")
    print(f"tracked: {rewrite_docx(tracked, True, wanted)} descriptions removed as {AUTHOR!r}")


if __name__ == "__main__":
    main()
