"""Minimal xlsx reader using only the standard library (no openpyxl in the rnaseq env).
read_sheets(path) -> {sheet_name: list of rows (list of str)}; first sheets only need shared strings and inline values."""
import zipfile, re
import xml.etree.ElementTree as ET
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}

def _col(ref):
    s = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in s:
        n = n * 26 + ord(ch) - 64
    return n - 1

def read_sheets(path, max_rows=None):
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", NS):
            ss.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    relmap = {r.get("Id"): r.get("Target") for r in rels}
    out = {}
    for sh in wb.find("m:sheets", NS):
        name = sh.get("name")
        target = relmap[sh.get("{%s}id" % NS["r"])]
        target = target.lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        root = ET.fromstring(z.read(target))
        rows = []
        for row in root.find("m:sheetData", NS):
            vals = {}
            for c in row.findall("m:c", NS):
                t = c.get("t"); v = c.find("m:v", NS)
                if t == "inlineStr":
                    txt = "".join(x.text or "" for x in c.iter("{%s}t" % NS["m"]))
                elif v is None:
                    txt = ""
                elif t == "s":
                    txt = ss[int(v.text)]
                else:
                    txt = v.text or ""
                vals[_col(c.get("r"))] = txt
            width = (max(vals) + 1) if vals else 0
            rows.append([vals.get(i, "") for i in range(width)])
            if max_rows and len(rows) >= max_rows:
                break
        out[name] = rows
    return out
