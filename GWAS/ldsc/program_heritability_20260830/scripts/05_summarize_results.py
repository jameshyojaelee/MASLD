#!/usr/bin/env python3
"""
Collect the partitioned heritability results and apply the pre-specified
multiple testing correction.

For every trait and every program the regression produces two numbers that mean
different things:

  Enrichment  the share of heritability the program's windows carry divided by
              the share of SNPs they contain. Easy to read, but unstable for
              small annotations and easily inflated when an annotation happens
              to overlap a baselineLD category such as coding sequence or
              conserved elements.

  Coefficient the per-SNP heritability the program adds on top of everything in
              the baselineLD model and on top of the cell-type background. This
              is the test that answers the question being asked. It is reported
              here in the standardised form tau*, which is the change in per-SNP
              heritability, as a fraction of total heritability, for a one
              standard deviation increase in the annotation
              (Gazal et al., Nature Genetics 2017):

                  tau* = tau * sd(annotation) * M / h2_g

              where M is the number of common reference SNPs the regression was
              normalised on and h2_g is the total heritability that same
              regression estimated.

Multiple testing: Benjamini-Hochberg across the 117 programs, separately within
each trait. The correction is applied to the one-sided coefficient p-value,
because only positive heritability enrichment is a meaningful finding. Traits in
different tiers are corrected separately and never pooled.

Usage: 05_summarize_results.py <results_root> <out_tsv> [ldscores_split_dir]
       results_root holds <trait>/<slug>.results and <trait>/<slug>.log
"""

import math
import os
import re
import sys

import numpy as np

REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE = os.path.join(REPO, "GWAS/ldsc/program_heritability_20260830")
MANIFEST = os.path.join(BASE, "gene_sets/annotation_manifest.tsv")
TRAITS = os.path.join(BASE, "config/traits.tsv")

H2_RE = re.compile(r"Total (?:Observed|Liability) scale h2: ([-\d.eE+]+) \(([-\d.eE+]+)\)")
SNP_RE = re.compile(r"After merging with regression SNP LD, (\d+) SNPs remain")

HEADLINE = (
    "No MASLD disease endpoint clears the heritability bar for stratified LD score\n"
    "regression. These results answer a question about serum ALT and AST, not about\n"
    "MASLD. A null on the disease endpoints would be uninformative, because those\n"
    "GWAS are too small to have shown anything either way. See config/traits.tsv.\n"
)


def normal_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def bh(pvalues):
    """Benjamini-Hochberg adjusted p-values, order preserved."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1.0)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def read_manifest():
    rows = {}
    with open(MANIFEST) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            rows[r["file_slug"]] = r
    return rows


def read_traits():
    rows = {}
    with open(TRAITS) as fh:
        header = None
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


def parse_results(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
    return rows


def annotation_snp_count(split_dir, slug):
    """Sum the annotation's common-SNP counts over the 22 autosomes."""
    total = 0.0
    for chrom in range(1, 23):
        path = os.path.join(split_dir, "%s.%d.l2.M_5_50" % (slug, chrom))
        with open(path) as fh:
            total += float(fh.read().split()[0])
    return total


def parse_h2(path):
    with open(path) as fh:
        text = fh.read()
    m = H2_RE.search(text)
    if not m:
        return None, None
    return float(m.group(1)), float(m.group(2))


def parse_regression_snps(path):
    with open(path) as fh:
        m = SNP_RE.search(fh.read())
    return int(m.group(1)) if m else None


def tier_caveat(tinfo):
    tier = tinfo["tier"]
    if tier == "primary":
        return ("liver injury biomarker, not a MASLD diagnosis; h2 z %s, attenuation %s"
                % (tinfo["h2_z"], tinfo["attenuation"]))
    if tier == "secondary":
        return ("h2 z %s is below the recommended bar of 7; treat as suggestive"
                % tinfo["h2_z"])
    if tier == "specificity_control":
        return "not a liver trait; included to test whether the programs are liver specific"
    return tinfo.get("reason", "")


def main():
    results_root, out_tsv = sys.argv[1], sys.argv[2]
    split_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(BASE, "ldscores_split")
    manifest = read_manifest()
    traits = read_traits()
    m_annot = {}

    records = []
    for trait, tinfo in traits.items():
        if tinfo["tier"] == "excluded":
            continue
        tdir = os.path.join(results_root, trait)
        if not os.path.isdir(tdir):
            continue
        for slug, minfo in manifest.items():
            if minfo["kind"] != "program":
                continue
            rpath = os.path.join(tdir, "%s.results" % slug)
            lpath = os.path.join(tdir, "%s.log" % slug)
            if not os.path.exists(rpath):
                continue
            rows = parse_results(rpath)
            # Look the program up by name. LDSC writes categories in the order
            # they were given to --ref-ld-chr and appends "L2_0" to each name,
            # but relying on position would fail silently if that ever changed.
            name = minfo["annotation"]
            hits = [r for r in rows if r["Category"].startswith(name)]
            if len(hits) != 1:
                sys.exit("%s: expected one category matching %r, found %d"
                         % (rpath, name, len(hits)))
            row = hits[0]
            trait_h2, trait_h2_se = parse_h2(lpath)
            prop_snps = float(row["Prop._SNPs"])
            tau = float(row["Coefficient"])
            tau_se = float(row["Coefficient_std_error"])
            z = float(row["Coefficient_z-score"])
            # LDSC normalises on M_5_50, the common reference SNPs. The
            # annotation's own count comes from its .l2.M_5_50 files, and the
            # total follows from the proportion LDSC reported.
            if slug not in m_annot:
                m_annot[slug] = annotation_snp_count(split_dir, slug)
            m_total = m_annot[slug] / prop_snps if prop_snps > 0 else float("nan")
            sd = math.sqrt(max(prop_snps * (1.0 - prop_snps), 0.0))
            tau_star = tau_star_se = float("nan")
            if trait_h2 and trait_h2 > 0 and m_total == m_total:
                scale = m_total * sd / trait_h2
                tau_star = tau * scale
                tau_star_se = tau_se * scale
            records.append({
                "trait": trait,
                "tier": tinfo["tier"],
                "tier_caveat": tier_caveat(tinfo),
                "n_regression_snps": parse_regression_snps(lpath),
                "program": minfo["annotation"],
                "lineage": minfo["lineage"],
                "n_genes_placed": minfo["n_genes_placed"],
                "frac_autosome": minfo["frac_autosome"],
                "prop_snps": prop_snps,
                "prop_h2": row["Prop._h2"],
                "prop_h2_se": row["Prop._h2_std_error"],
                "enrichment": row["Enrichment"],
                "enrichment_se": row["Enrichment_std_error"],
                "enrichment_p": row["Enrichment_p"],
                "tau": tau,
                "tau_se": tau_se,
                "tau_z": z,
                "tau_p_one_sided": normal_sf(z),
                "tau_star": tau_star,
                "tau_star_se": tau_star_se,
                "model_h2": trait_h2,
                "model_h2_se": trait_h2_se,
                "m_annotation_5_50": m_annot[slug],
                "m_total_5_50": m_total,
            })

    if not records:
        sys.exit("no .results files found under %s" % results_root)

    # Benjamini-Hochberg within trait, across the 117 programs.
    by_trait = {}
    for i, r in enumerate(records):
        by_trait.setdefault(r["trait"], []).append(i)
    for trait, idx in by_trait.items():
        adj = bh([records[i]["tau_p_one_sided"] for i in idx])
        adj_e = bh([float(records[i]["enrichment_p"]) for i in idx])
        for i, a, ae in zip(idx, adj, adj_e):
            records[i]["tau_p_bh_within_trait"] = a
            records[i]["enrichment_p_bh_within_trait"] = ae
            records[i]["n_tests_in_family"] = len(idx)

    cols = ["trait", "tier", "tier_caveat", "program", "lineage",
            "n_genes_placed", "frac_autosome", "n_regression_snps",
            "prop_snps", "prop_h2", "prop_h2_se", "enrichment", "enrichment_se",
            "enrichment_p", "enrichment_p_bh_within_trait", "tau", "tau_se", "tau_z",
            "tau_p_one_sided", "tau_p_bh_within_trait", "tau_star", "tau_star_se",
            "model_h2", "model_h2_se", "m_annotation_5_50", "m_total_5_50",
            "n_tests_in_family"]
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(records, key=lambda r: (r["tier"], r["trait"], r["tau_p_one_sided"])):
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    outdir = os.path.dirname(os.path.abspath(out_tsv))
    with open(os.path.join(outdir, "HEADLINE.txt"), "w") as fh:
        fh.write(HEADLINE)
    with open(os.path.join(outdir, "excluded_traits.tsv"), "w") as fh:
        fh.write("trait\tdescription\th2_z\tattenuation\texclusion_reason\n")
        for t, ti in sorted(traits.items()):
            if ti["tier"] == "excluded":
                fh.write("\t".join([t, ti["description"], ti["h2_z"],
                                     ti["attenuation"], ti["reason"]]) + "\n")

    print(HEADLINE)
    print("wrote %s (%d rows)" % (out_tsv, len(records)))
    for trait, idx in sorted(by_trait.items()):
        n = sum(1 for i in idx if records[i]["tau_p_bh_within_trait"] < 0.05)
        ne = sum(1 for i in idx if records[i]["enrichment_p_bh_within_trait"] < 0.05)
        print("  %-16s %3d programs tested, %2d pass BH on the coefficient, "
              "%2d pass BH on enrichment" % (trait, len(idx), n, ne))


if __name__ == "__main__":
    main()
