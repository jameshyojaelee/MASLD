#!/usr/bin/env python3
"""
Report the cell-type background annotations on their own, before any program
result is read.

Each background is the set of every gene that entered the Hotspot
autocorrelation test in one cell type, put through the same 100 kb window rule
as the programs. Together they cover 50 to 68 percent of the autosome, and every
program's windows sit inside its own background.

If a background is itself strongly enriched for a trait, then a program inside
it starts from an elevated baseline, and its own enrichment has to be read
against that rather than against one. This is worth knowing before the program
results are opened, not after, which is why it runs first and is reported
separately.

Usage: 06_summarize_backgrounds.py <results_root> <out_tsv>
"""

import math
import os
import re
import sys

REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE = os.path.join(REPO, "GWAS/ldsc/program_heritability_20260830")
MANIFEST = os.path.join(BASE, "gene_sets/annotation_manifest.tsv")
TRAITS = os.path.join(BASE, "config/traits.tsv")

H2_RE = re.compile(r"Total (?:Observed|Liability) scale h2: ([-\d.eE+]+) \(([-\d.eE+]+)\)")
SNP_RE = re.compile(r"After merging with regression SNP LD, (\d+) SNPs remain")

LINEAGES = ["hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells"]


def read_traits():
    rows, header = {}, None
    with open(TRAITS) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            r = dict(zip(header, fields))
            rows[r["trait"]] = r
    return rows


def read_manifest():
    rows, header = {}, None
    with open(MANIFEST) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            rows[r["file_slug"]] = r
    return rows


def parse_log(path):
    with open(path) as fh:
        text = fh.read()
    h2 = H2_RE.search(text)
    snps = SNP_RE.search(text)
    return (float(h2.group(1)) if h2 else None,
            float(h2.group(2)) if h2 else None,
            int(snps.group(1)) if snps else None)


def main():
    results_root, out_tsv = sys.argv[1:3]
    traits = read_traits()
    manifest = read_manifest()

    cols = ["trait", "tier", "tier_caveat", "background", "lineage",
            "n_genes_listed", "n_genes_placed", "frac_autosome",
            "prop_snps", "prop_h2", "prop_h2_se", "enrichment", "enrichment_se",
            "enrichment_p", "tau", "tau_se", "tau_z", "tau_p_one_sided",
            "model_h2", "model_h2_se", "n_regression_snps"]
    out = []

    for trait, tinfo in sorted(traits.items()):
        if tinfo["tier"] == "excluded":
            continue
        tdir = os.path.join(results_root, trait)
        for lineage in LINEAGES:
            slug = "background-%s" % lineage
            rpath = os.path.join(tdir, "_background_only_%s.results" % lineage)
            lpath = os.path.join(tdir, "_background_only_%s.log" % lineage)
            if not os.path.exists(rpath):
                continue
            with open(rpath) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                rows = [dict(zip(header, l.rstrip("\n").split("\t")))
                        for l in fh if l.strip()]
            name = manifest[slug]["annotation"]
            hits = [r for r in rows if r["Category"].startswith(name)]
            if len(hits) != 1:
                sys.exit("%s: expected one category matching %r, found %d"
                         % (rpath, name, len(hits)))
            row = hits[0]
            h2, h2se, nsnp = parse_log(lpath)
            z = float(row["Coefficient_z-score"])
            out.append({
                "trait": trait, "tier": tinfo["tier"],
                "tier_caveat": tier_caveat(tinfo["tier"], tinfo),
                "background": name, "lineage": lineage,
                "n_genes_listed": manifest[slug]["n_genes_listed"],
                "n_genes_placed": manifest[slug]["n_genes_placed"],
                "frac_autosome": manifest[slug]["frac_autosome"],
                "prop_snps": row["Prop._SNPs"], "prop_h2": row["Prop._h2"],
                "prop_h2_se": row["Prop._h2_std_error"],
                "enrichment": row["Enrichment"],
                "enrichment_se": row["Enrichment_std_error"],
                "enrichment_p": row["Enrichment_p"],
                "tau": row["Coefficient"], "tau_se": row["Coefficient_std_error"],
                "tau_z": z, "tau_p_one_sided": 0.5 * math.erfc(z / math.sqrt(2.0)),
                "model_h2": h2, "model_h2_se": h2se, "n_regression_snps": nsnp,
            })

    if not out:
        sys.exit("no background results found under %s" % results_root)

    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    print("wrote %s (%d rows)" % (out_tsv, len(out)))
    print()
    print("Background annotations on their own, against baselineLD v2.2.")
    print("Read every program enrichment for a cell type against its background here.")
    print("%-16s %-16s %8s %10s %10s %8s" % ("trait", "lineage", "propSNP", "enrichment", "enrich_p", "tau_z"))
    for r in out:
        print("%-16s %-16s %8.4f %10.4f %10s %8.2f"
              % (r["trait"], r["lineage"], float(r["prop_snps"]),
                 float(r["enrichment"]), r["enrichment_p"], r["tau_z"]))


def tier_caveat(tier, tinfo):
    if tier == "primary":
        return ("liver injury biomarker, not a MASLD diagnosis; "
                "h2 z %s, attenuation %s" % (tinfo["h2_z"], tinfo["attenuation"]))
    if tier == "secondary":
        return ("h2 z %s is below the recommended bar of 7; treat as suggestive"
                % tinfo["h2_z"])
    if tier == "specificity_control":
        return "not a liver trait; included to test whether the programs are liver specific"
    return tinfo.get("reason", "")


if __name__ == "__main__":
    main()
