#!/usr/bin/env python3
"""
69a_borzoi_magnitude_prep.py  --  Stage-2b(Axis-4) LEAD SELECTION for the Borzoi
MAGNITUDE / gene-attribution nominator (direction-INDEPENDENT).

ADDITIVE ONLY -- new file. Does NOT edit 46d/78/27a or any existing 60-65 script.
Consumes only finemapping outputs; emits a lead table for 69_borzoi_magnitude_nom.py.

WHY THIS STAGE EXISTS
  The eQTL-DIRECTION gate FAILED for Borzoi AND AlphaGenome (auROC ~0.53-0.56), so
  direction is off the table. But Borzoi |logSED| tracks |eQTL beta| (rho 0.39,
  significant) -- a legitimate MAGNITUDE nomination (the published Sniff use). This
  stage picks, for each independent eQTL-ABSENT fine-mapped signal, the single
  strongest-PIP regulatory SNV to hand to Borzoi as the lead. Gene attribution
  (enumerate candidates in-window, score each, report best-magnitude gene) happens
  in 69_borzoi_magnitude_nom.py.

TARGET SET (eQTL-absent): substrate rows with colocalizes==FALSE & var_class==
  regulatory & is_indel==FALSE & a resolved hg38 SNV (ref_hg38 non-empty). These are
  the ~84% of fine-mapped liver signals whose effector gene does NOT colocalize.

SIGNAL / LOCUS DEFINITION (reproducible, GWAS-complete)
  merged_loci_map.csv only harmonizes a subset of the 55-GWAS portfolio (241 of the
  study-loci; ~18k eQTL-absent credible-set rows fall in studies it never mapped), so
  collapsing on its merge_group would silently drop most signals. Instead we
  distance-clump the eQTL-absent regulatory leads themselves on hg38 coordinates:

    1. For every eQTL-absent regulatory SNV, take its STRONGEST fine-mapping PIP
       anywhere in the 55-GWAS credible_sets (max recommended_pip over all study rows
       for that variant_id; recommended_pip is the SuSiE/CARMA-reconciled PIP).
    2. Keep variants that are genuinely fine-mapped: lead PIP >= PIP_FLOOR (0.10 =
       in a 95% credible set with non-trivial mass). max_pip==0 consequence-only rows
       are dropped.
    3. Single-linkage clump the survivors within CLUMP_BP (250 kb ~ Borzoi half
       receptive field) on the same hg38 contig; each clump = one physical signal.
    4. Signal lead = the max-PIP member (ties -> higher hg38 pos, deterministic).

  250 kb was chosen so that within one signal the per-candidate-gene Borzoi windows
  substantially overlap (two variants >250 kb apart see largely disjoint gene sets and
  are treated as separate signals). This is LD-agnostic distance clumping on already
  fine-mapped credible-set variants, not a re-finemapping.

OUTPUT: results/seqfunc/borzoi_magnitude_leads.tsv
  cols: locus_id, variant_id_hg19, chr, pos_hg38, ref_hg38, alt_hg38, lead_pip,
        gene_substrate, ensembl_substrate, n_cs_variants_in_clump, clump_span_bp
"""
import csv
import os
import sys

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM = os.path.join(BASE, "GWAS/finemapping/results")
SUBSTRATE = os.path.join(FM, "seqfunc/variant_substrate_hg38.tsv")
CREDSETS = os.path.join(FM, "credible_sets.csv")
OUT = os.path.join(FM, "seqfunc/borzoi_magnitude_leads.tsv")

PIP_FLOOR = float(os.environ.get("MAG_PIP_FLOOR", "0.10"))
CLUMP_BP = int(os.environ.get("MAG_CLUMP_BP", "250000"))


def log(m):
    print(m, flush=True)


def main():
    # --- 1. eQTL-absent regulatory SNV substrate (hg38-resolved) -------------
    sub = {}
    with open(SUBSTRATE) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            if (x["colocalizes"] == "FALSE" and x["var_class"] == "regulatory"
                    and x["is_indel"] == "FALSE" and x["ref_hg38"]):
                sub[x["variant_id_hg19"]] = x
    log(f"[1] eQTL-absent regulatory hg38-SNV substrate variants: {len(sub)}")

    # --- 2. strongest fine-mapping PIP per variant across the 55-GWAS credsets
    best_pip = {}
    with open(CREDSETS) as f:
        for x in csv.DictReader(f):
            vid = x["variant_id"]
            if vid not in sub:
                continue
            try:
                pip = float(x["recommended_pip"])
            except (ValueError, TypeError):
                try:
                    pip = float(x["max_pip"])
                except (ValueError, TypeError):
                    pip = 0.0
            if vid not in best_pip or pip > best_pip[vid]:
                best_pip[vid] = pip
    log(f"[2] eQTL-absent variants found in credible_sets: {len(best_pip)}")

    # --- 3. fine-mapped survivors (lead PIP >= floor) ------------------------
    leads = []
    for vid, s in sub.items():
        pip = best_pip.get(vid, 0.0)
        if pip < PIP_FLOOR:
            continue
        leads.append({
            "variant_id_hg19": vid,
            "chr": s["chr"],
            "pos_hg38": int(s["pos_hg38"]),
            "ref_hg38": s["ref_hg38"],
            "alt_hg38": s["alt_hg38"],
            "pip": pip,
            "gene_substrate": s["gene"],
            "ensembl_substrate": s["ensembl"],
        })
    log(f"[3] fine-mapped eQTL-absent reg SNVs (PIP>={PIP_FLOOR}): {len(leads)}")
    if not leads:
        log("FATAL: no leads survived the PIP floor.")
        sys.exit(2)

    # --- 4. single-linkage distance clump per contig -------------------------
    by_chr = {}
    for L in leads:
        by_chr.setdefault(L["chr"], []).append(L)

    signals = []
    for chrom, rows in by_chr.items():
        rows.sort(key=lambda r: r["pos_hg38"])
        cur = [rows[0]]
        for r in rows[1:]:
            if r["pos_hg38"] - cur[-1]["pos_hg38"] <= CLUMP_BP:
                cur.append(r)
            else:
                signals.append(cur)
                cur = [r]
        signals.append(cur)

    # --- 5. per-signal lead = max PIP (tie -> higher pos) --------------------
    out_rows = []
    for clump in signals:
        lead = max(clump, key=lambda r: (r["pip"], r["pos_hg38"]))
        span = max(r["pos_hg38"] for r in clump) - min(r["pos_hg38"] for r in clump)
        locus_id = f"{lead['chr']}:{lead['pos_hg38']}"
        out_rows.append({
            "locus_id": locus_id,
            "variant_id_hg19": lead["variant_id_hg19"],
            "chr": lead["chr"],
            "pos_hg38": lead["pos_hg38"],
            "ref_hg38": lead["ref_hg38"],
            "alt_hg38": lead["alt_hg38"],
            "lead_pip": round(lead["pip"], 6),
            "gene_substrate": lead["gene_substrate"],
            "ensembl_substrate": lead["ensembl_substrate"],
            "n_cs_variants_in_clump": len(clump),
            "clump_span_bp": span,
        })
    out_rows.sort(key=lambda r: (r["chr"], r["pos_hg38"]))

    cols = ["locus_id", "variant_id_hg19", "chr", "pos_hg38", "ref_hg38", "alt_hg38",
            "lead_pip", "gene_substrate", "ensembl_substrate",
            "n_cs_variants_in_clump", "clump_span_bp"]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    pips = sorted((r["lead_pip"] for r in out_rows), reverse=True)
    log(f"[4-5] distinct eQTL-absent signals (clump {CLUMP_BP//1000}kb): {len(out_rows)}")
    log(f"      lead PIP  >=0.9:{sum(p>=0.9 for p in pips)}  "
        f">=0.5:{sum(p>=0.5 for p in pips)}  >=0.1:{sum(p>=0.1 for p in pips)}")
    log(f"      wrote {OUT} ({len(out_rows)} signals)")


if __name__ == "__main__":
    main()
