#!/usr/bin/env python3
"""N5 gate: which perturbations of HSD17B13 or PNPLA3 exist in the data on disk.

Reads metadata only (signature annotation tables, compound annotation, the Saunders
target list, the H3 per-gene tables). No expression value of any perturbation is read.
Run with: python -I n5_gate_presence_check.py <out_dir>
"""
import csv
import gzip
import hashlib
import os
import sys
from collections import Counter

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
LINCS = f"{ROOT}/Analysis/Hypotheses/H47-20261008T214453Z/downloads_untrusted"
SAUNDERS = f"{ROOT}/Analysis/Hypotheses/H48-20261008T014825Z/step1/knockout_eligibility.tsv"
H3 = f"{ROOT}/Analysis/Hypotheses/H3-20261009T180331Z/results"
ALL_GENES = {
    "PNPLA3_G_dosage": f"{H3}/signature-pnpla3-20261009T184935Z/all_genes_PNPLA3_G_dosage.tsv.gz",
    "HSD17B13_TA_carrier": f"{H3}/signature-hsd-20261009T184935Z/all_genes_HSD17B13_TA_carrier.tsv.gz",
}
SIG_INFO = {
    "GSE92742": f"{LINCS}/GSE92742_Broad_LINCS_sig_info.txt.gz",
    "GSE70138": f"{LINCS}/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
}
COMPOUND_INFO = f"{LINCS}/compoundinfo_beta.txt"
GENE_INFO = f"{LINCS}/GSE92742_Broad_LINCS_gene_info.txt.gz"
LIVER_LINES = {"HEPG2", "HUH7", "PHH"}
# The two genes of interest first; the rest are family members listed as context.
FOCUS = ["HSD17B13", "PNPLA3"]
CONTEXT = ["HSD17B11", "HSD17B12", "HSD17B4", "PNPLA2", "PNPLA5", "PNPLA8"]
SAUNDERS_FOCUS = {"Hsd17b13", "Pnpla3", "Pnpla2", "Mtarc1", "Gpam", "Cideb", "Tm6sf2", "Gckr", "Mboat7"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def main(out_dir):
    os.makedirs(out_dir, exist_ok=False)
    with open(os.path.join(out_dir, "input_sha256.tsv"), "w") as fh:
        fh.write("path\tsha256\n")
        for p in [*SIG_INFO.values(), COMPOUND_INFO, GENE_INFO, SAUNDERS, *ALL_GENES.values()]:
            fh.write(f"{p}\t{sha256(p)}\n")

    # 1. LINCS signature records by gene, perturbation type and cell line
    rows = []
    liver_pert = Counter()
    liver_genes = {}
    for series, path in SIG_INFO.items():
        with gzip.open(path, "rt") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for r in rdr:
                name, ptype, cell = r["pert_iname"], r["pert_type"], r["cell_id"]
                if name in FOCUS or name in CONTEXT:
                    rows.append((series, name, ptype, cell, 1))
                if cell in LIVER_LINES:
                    liver_pert[(series, ptype)] += 1
                    if ptype.startswith("trt_sh") or ptype in ("trt_oe", "trt_xpr"):
                        liver_genes.setdefault((series, ptype), set()).add(name)
    agg = Counter()
    for series, name, ptype, cell, n in rows:
        agg[(series, name, ptype, cell)] += n
    with open(os.path.join(out_dir, "lincs_records_focus_genes.tsv"), "w") as fh:
        fh.write("series\tgene\tpert_type\tcell_id\tn_signatures\tliver_line\tgene_class\n")
        for (series, name, ptype, cell), n in sorted(agg.items()):
            cls = "focus" if name in FOCUS else "context"
            fh.write(f"{series}\t{name}\t{ptype}\t{cell}\t{n}\t{cell in LIVER_LINES}\t{cls}\n")
    with open(os.path.join(out_dir, "lincs_focus_summary.tsv"), "w") as fh:
        fh.write("gene\tn_signatures_any_series_any_line\tn_signatures_liver_lines\n")
        for g in FOCUS + CONTEXT:
            tot = sum(n for (s, name, pt, c), n in agg.items() if name == g)
            liv = sum(n for (s, name, pt, c), n in agg.items() if name == g and c in LIVER_LINES)
            fh.write(f"{g}\t{tot}\t{liv}\n")
    with open(os.path.join(out_dir, "lincs_liver_line_perturbation_context.tsv"), "w") as fh:
        fh.write("series\tpert_type\tn_signatures_liver_lines\tn_distinct_perturbed_genes_liver_lines\n")
        for (series, ptype), n in sorted(liver_pert.items()):
            ng = len(liver_genes.get((series, ptype), ()))
            fh.write(f"{series}\t{ptype}\t{n}\t{ng if ng else ''}\n")

    # 2. compound annotation: any compound whose target field names a focus gene
    hits = []
    n_compounds = 0
    with open(COMPOUND_INFO) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            n_compounds += 1
            tokens = {t.strip().upper() for t in (r.get("target") or "").replace(",", "|").split("|")}
            aliases = (r.get("compound_aliases") or "").upper()
            for g in FOCUS:
                if g in tokens or g in aliases:
                    hits.append((g, r["cmap_name"], r.get("target", ""), r.get("moa", "")))
    with open(os.path.join(out_dir, "compound_annotation_focus_targets.tsv"), "w") as fh:
        fh.write("gene\tcmap_name\ttarget\tmoa\n")
        for h in hits:
            fh.write("\t".join(h) + "\n")
        fh.write(f"# compounds_in_table\t{n_compounds}\t\t\n")

    # 3. Saunders library
    with open(SAUNDERS) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        targets = [r for r in rdr]
    with open(os.path.join(out_dir, "saunders_focus_targets.tsv"), "w") as fh:
        fh.write("target_gene\tin_library\ttested\n")
        present = {r["target_gene"]: r["tested"] for r in targets}
        for g in sorted(SAUNDERS_FOCUS):
            fh.write(f"{g}\t{g in present}\t{present.get(g, '')}\n")
        fh.write(f"# targets_in_library\t{len(targets)}\t{sum(r['tested'] == 'True' for r in targets)}\n")

    # 4. signature composition against LINCS landmarks
    landmarks = set()
    with gzip.open(GENE_INFO, "rt") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            if r["pr_is_lm"] == "1":
                landmarks.add(r["pr_gene_symbol"])
    with open(os.path.join(out_dir, "signature_landmark_composition.tsv"), "w") as fh:
        fh.write("model\tgenes_tested\tlandmark_genes_in_table\tlandmark_genes_fdr05\tlandmark_up_t_gt0\tlandmark_down_t_lt0\tn_landmarks_lincs\n")
        for model, path in ALL_GENES.items():
            n = lm = lm_fdr = up = down = 0
            with gzip.open(path, "rt") as gh:
                rdr = csv.DictReader(gh, delimiter="\t")
                for r in rdr:
                    n += 1
                    if r["symbol"] in landmarks:
                        lm += 1
                        t = float(r["t"])
                        up += t > 0
                        down += t < 0
                        lm_fdr += float(r["adj.P.Val"]) < 0.05
            fh.write(f"{model}\t{n}\t{lm}\t{lm_fdr}\t{up}\t{down}\t{len(landmarks)}\n")
    print("done", out_dir)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python -I n5_gate_presence_check.py <out_dir>")
    main(sys.argv[1])
