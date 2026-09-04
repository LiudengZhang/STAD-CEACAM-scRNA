import sys, zipfile, re
from xml.etree import ElementTree as ET
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
def paras(path):
    z=zipfile.ZipFile(path)
    root=ET.fromstring(z.read('word/document.xml'))
    out=[]
    for p in root.iter('{%s}p'%NS['w']):
        t=''.join(n.text or '' for n in p.iter('{%s}t'%NS['w']))
        out.append(t)
    return out
if __name__=='__main__':
    for i,t in enumerate(paras(sys.argv[1]),1):
        if len(sys.argv)>2:
            if not re.search(sys.argv[2], t, re.I): continue
        print(f"[{i}] {t}")
