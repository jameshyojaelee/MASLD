"""Independently compare rendered numbers and sources with adopted per-gene calls."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

import pandas as pd
import pdfplumber
from pdfminer.pdftypes import resolve1

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = ROOT / "Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures"
CAT = ROOT / "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106"
index = pd.read_parquet(CAT / "passport_gene_index.parquet")
e = pd.read_parquet(CAT / "passport_evidence_long.parquet")
g = e[e.evidence_domain == "genetics"].set_index("ensembl_id").call_state
r = e[e.evidence_domain == "transcriptomics"].set_index("ensembl_id").call_state
matrix = pd.crosstab(g, r)
src = pd.read_csv(OUT / "source_tables/6B.tsv", sep="\t")
for row in src[src.record == "cell"].itertuples():
    assert matrix.loc[row.genetics_call, row.rna_call] == row.n_genes
assert int(matrix.to_numpy().sum()) == 28233
assert int(matrix.loc["untestable", "untestable"]) == 8
routes = pd.read_csv(OUT / "source_tables/6D.tsv", sep="\t")
bars = routes[(routes.hatched == False) & (routes.bar_left.notna())]
classes = bars.groupby("primary_evidence_class").n_genes.sum().astype(int)
assert classes.to_dict() == index.primary_evidence_class.value_counts().to_dict()
native = pd.read_csv(OUT / "source_tables/fig3g.tsv", sep="\t")
saved = pd.read_csv(ROOT / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/gene_level_20261021T165252Z/h57_universe_rule_a_sens.tsv", sep="\t")
assert set(native.g) == set(saved.g) and len(native) == 3589
assert native.display_class.value_counts().to_dict() == {"Neither": 2980, "Genetics": 400, "Transcriptomics": 183, "Both": 26}
summary = pd.read_csv(OUT / "source_tables/native_overlap_summary.tsv", sep="\t").iloc[0]
assert (int(summary.N), int(summary.native_interface), int(summary.catalog_jointly_eligible)) == (3589, 14704, 14502)
records = []
for path in sorted((OUT / "panels").glob("*.pdf")):
    text = subprocess.check_output(["pdftotext", "-layout", str(path), "-"], text=True)
    (OUT / "proofs" / (path.stem + ".txt")).write_text(text)
    fonts = []
    with pdfplumber.open(path) as font_pdf:
        for name, ref in resolve1(font_pdf.pages[0].page_obj.resources["Font"]).items():
            font = resolve1(ref)
            descendants = resolve1(font.get("DescendantFonts", []))
            actual = [resolve1(x) for x in descendants] or [font]
            for item in actual:
                assert "Type3" not in str(item.get("Subtype")), (path.name, item)
                descriptor = resolve1(item.get("FontDescriptor", {}))
                assert any(k in descriptor for k in ("FontFile", "FontFile2", "FontFile3")), (path.name, name, "unembedded font")
                fonts.append((name, str(item.get("Subtype")), str(item.get("BaseFont"))))
    (OUT / "proofs" / (path.stem + "_fonts.txt")).write_text(repr(fonts) + "\n")
    with pdfplumber.open(path) as pdf:
        assert len(pdf.pages) == 1
        page = pdf.pages[0]
        sizes = sorted({round(c["size"], 2) for c in page.chars if c["upright"] and c["text"].strip()})
        assert sizes == [6.0], (path.name, sizes)
        overflow = [c["text"] for c in page.chars if c["x0"] < -.15 or c["x1"] > page.width + .15 or c["top"] < -.15 or c["bottom"] > page.height + .15]
        assert not overflow, (path.name, "text outside page", overflow)
        assert not any("Bold" in c["fontname"] for c in page.chars)
        records.append(dict(panel=path.name, width_pt=page.width, height_pt=page.height,
                            text_sizes=sizes, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    if path.name.startswith("fig6b"):
        for value in ["12,958", "1,118", "9,507", "4,214", "8,790", "709", "195"]:
            assert value in text, value
        assert "none by construction" not in text
    if path.name.startswith("fig6d"):
        for value in ["13,731", "12,958", "1,118", "9,499", "719"]:
            assert value in text, value
    if path.name.startswith("fig3g"):
        assert "2,980" in text and "183" in text and "3,149" not in text
    if path.name.startswith("fig1b"):
        assert "3,589" in text and "209" in text and "21.4" in text
        assert "Expression-matched non-colocalized" not in text
    subprocess.run(["pdftoppm", "-scale-to", "1600", "-png", "-singlefile", str(path), str(OUT / "proofs" / path.stem)], check=True)
(OUT / "checks.json").write_text(json.dumps(dict(passed=True, per_gene_catalog_matrix=True,
    per_gene_routes=True, native_ids=True, font_and_geometry_checks=True, panels=records), indent=2) + "\n")
pd.DataFrame(records).to_csv(OUT / "panel_manifest.tsv", sep="\t", index=False)
print("PASS: adopted per-gene counts/routes, native eligibility, rendered values, 6pt text, page bounds and no Type3 fonts")
