"""Sentence-level extractor for NF-kB / TNFa-signalling claims.

Read-only. Splits a text source into sentences and keeps those that match the
NF-kB claim vocabulary. Used for CLAIM_COMPARISON.md and for its positive control.
"""
import re, sys, zipfile
from xml.etree import ElementTree as ET

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'

# The claim vocabulary. A sentence is an NF-kB / TNFa-signalling claim candidate
# if it names the pathway, the Hallmark set, or the GSEA statistic attached to it.
PAT = re.compile(
    r"NF\s*[-‐-―−]?\s*(?:κB|kB|KB|ΚB)"     # NF-kB / NF-κB
    r"|NFKB[12]?|NFKB\d"
    r"|TNF\s*[αa]?\s*(?:signaling|signalling)"
    r"|TNFα|TNF-α|TNF-alpha"
    r"|Hallmark"
    r"|top[- ]ranked"
    r"|NES\s*=",
    re.I)

SPLIT = re.compile(r'(?<=[.;:!?])\s+(?=[A-Z0-9(“"“])')

def docx_paras(path):
    z = zipfile.ZipFile(path)
    root = ET.fromstring(z.read('word/document.xml'))
    return [''.join(n.text or '' for n in p.iter(f'{{{W}}}t'))
            for p in root.iter(f'{{{W}}}p')]

def sentences(para):
    para = para.strip()
    if not para:
        return []
    return [s.strip() for s in SPLIT.split(para) if s.strip()]

def scan_paras(paras, label):
    hits = []
    for i, p in enumerate(paras, 1):
        for j, s in enumerate(sentences(p), 1):
            if PAT.search(s):
                hits.append((label, i, j, s))
    return hits

def main():
    src = sys.argv[1]
    label = sys.argv[2] if len(sys.argv) > 2 else src
    if src.endswith('.docx'):
        paras = docx_paras(src)
    else:
        paras = open(src, encoding='utf-8').read().split('\n')
    for h in scan_paras(paras, label):
        print(f"{h[0]}\tpara{h[1]}\ts{h[2]}\t{h[3]}")

if __name__ == '__main__':
    main()
