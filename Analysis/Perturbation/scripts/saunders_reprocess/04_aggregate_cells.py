#!/usr/bin/env python3
"""
Aggregate per-SRR counts into per-cell sgRNA assignments.

Steps:
1. Read all <SRR>.counts.tsv.gz files
2. Aggregate (cell_barcode_16, sample_bc) → {sgrna_id: (n_reads, n_umis)}
   (sublib identity is inferred from which SRR the reads came from; multiple SRRs
    can map to the same (CB, sample_bc) if a cell is observed in different SRRs)
3. For each cell × sample_bc:
     - Compute dominant_fraction = top_sgrna_umis / sum(umis)
     - Bootstrap test (pftools-style): is top sgRNA significantly enriched vs random?
4. Emit per-cell assignment:
     cell_barcode_24 = CB16 + sample_bc8 (matches h5 format)
     fields: assigned sgrna_id, target gene, CR, dominant_fraction, n_umis_top,
             n_umis_total, padj_bootstrap, n_sgrnas_observed, confidence_tier
5. Optional: intersect with the GEO h5 cell barcodes
"""
import argparse, csv, glob, gzip, json, os
from collections import defaultdict, Counter
import numpy as np
from scipy.stats import binomtest


def load_all_counts(in_dir):
    """Returns dict[(cb16, sample_bc)][sgrna_id] = (n_reads, n_umis)."""
    cells = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    gene_for_sg = {}
    cr_for_sg = {}
    files = sorted(glob.glob(os.path.join(in_dir, "*.counts.tsv.gz")))
    print(f"Found {len(files)} per-SRR count files")
    for fp in files:
        with gzip.open(fp, "rt") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for row in rd:
                cb = row["cell_barcode"]
                pbc = row["sample_bc"]
                sgid = row["sgrna_id"]
                n_reads = int(row["n_reads"])
                n_umis = int(row["n_umis"])
                cells[(cb, pbc)][sgid][0] += n_reads
                cells[(cb, pbc)][sgid][1] += n_umis
                gene_for_sg[sgid] = row["gene"]
                cr_for_sg[sgid] = row["CR"]
    return cells, gene_for_sg, cr_for_sg


def assign_per_cell(cell_sg_counts, min_umis=3, min_dominant_frac=0.5, min_top_umis=2):
    """For each cell, decide assigned sgRNA and confidence tier."""
    if not cell_sg_counts:
        return None
    sgs = list(cell_sg_counts.items())  # [(sgid, [n_reads, n_umis]), ...]
    sgs.sort(key=lambda x: -x[1][1])    # sort by n_umis desc
    total_umis = sum(s[1][1] for s in sgs)
    total_reads = sum(s[1][0] for s in sgs)
    if total_umis < min_umis:
        return {"tier": "low_total_umis", "n_sgrnas": len(sgs), "total_umis": total_umis,
                "total_reads": total_reads}
    top_sgid, (top_reads, top_umis) = sgs[0]
    dominant_frac = top_umis / total_umis
    # binomial test: is top sgRNA's UMI fraction > 1/N_sgrnas (uniform null)?
    # Use the simplest possible: enrichment vs uniform expectation
    # n_sgrnas in library = 456; null = 1/456 = 0.0022
    p = binomtest(top_umis, total_umis, p=1.0/456, alternative="greater").pvalue
    tier = "confident"
    if dominant_frac < min_dominant_frac:
        tier = "ambiguous_low_dominance"
    elif top_umis < min_top_umis:
        tier = "low_top_umi"
    return {
        "tier": tier,
        "top_sgid": top_sgid,
        "top_umis": top_umis,
        "top_reads": top_reads,
        "total_umis": total_umis,
        "total_reads": total_reads,
        "dominant_frac": dominant_frac,
        "n_sgrnas": len(sgs),
        "second_sgid": sgs[1][0] if len(sgs) > 1 else None,
        "second_umis": sgs[1][1][1] if len(sgs) > 1 else 0,
        "binom_p": p,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--h5-glob",
                    default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/datasets/saunders2025/raw/h5/GSM8478327_*_filtered_feature_bc_matrix.h5",
                    help="If provided, intersect assignments with these h5 cell barcodes.")
    ap.add_argument("--min-umis", type=int, default=3)
    ap.add_argument("--min-dominant-frac", type=float, default=0.5)
    args = ap.parse_args()

    print("Loading per-SRR counts...")
    cells, gene_for_sg, cr_for_sg = load_all_counts(args.in_dir)
    print(f"Total unique (CB, sample_bc) cells: {len(cells):,}")

    # Per-cell assignment
    print("Assigning sgRNAs per cell...")
    assigns = {}
    tier_counter = Counter()
    for key, sg_counts in cells.items():
        a = assign_per_cell(sg_counts, args.min_umis, args.min_dominant_frac)
        if a is not None:
            assigns[key] = a
            tier_counter[a["tier"]] += 1

    print(f"Tier distribution: {dict(tier_counter)}")

    # Optional h5 intersection
    h5_barcodes = set()
    if args.h5_glob:
        try:
            import h5py
            h5_files = sorted(glob.glob(args.h5_glob))
            print(f"Loading barcodes from {len(h5_files)} h5 files for intersection")
            for hp in h5_files:
                with h5py.File(hp, "r") as h:
                    for bc in h["matrix"]["barcodes"][:]:
                        b = bc.decode()
                        # bc format: <16-nt CB><8-nt sample_bc>-1
                        b_strip = b.rstrip("-1").split("-")[0]
                        h5_barcodes.add(b_strip)
            print(f"Total h5 cells: {len(h5_barcodes):,}")
        except Exception as e:
            print(f"WARN: failed to load h5 barcodes: {e}")

    # Map BC001..BC004 -> raw 8-mer (h5 cell barcode uses the raw 8-mer)
    BC_TO_RAW = {"BC001": "ACTTTAGG", "BC002": "AACGGGAA",
                 "BC003": "AGTAGGCT", "BC004": "ATGTTGAC"}

    # Write output
    print(f"Writing {args.out}")
    with gzip.open(args.out, "wt") as out:
        out.write("\t".join([
            "cell_barcode_24", "cell_barcode_16", "sample_bc",
            "tier", "in_h5",
            "sgrna_id", "target_gene", "CR",
            "n_umis", "n_reads", "total_umis", "total_reads",
            "dominant_frac", "n_sgrnas", "second_sgid", "second_umis", "binom_p",
        ]) + "\n")
        n_total = 0
        n_in_h5 = 0
        n_confident_in_h5 = 0
        for (cb, pbc), a in assigns.items():
            raw_bc = BC_TO_RAW.get(pbc, pbc)
            cb24 = cb + raw_bc
            in_h5 = "TRUE" if cb24 in h5_barcodes else "FALSE"
            if cb24 in h5_barcodes:
                n_in_h5 += 1
                if a["tier"] == "confident":
                    n_confident_in_h5 += 1
            n_total += 1
            sg = a.get("top_sgid", "")
            out.write("\t".join([
                cb24, cb, pbc,
                a["tier"], in_h5,
                sg, gene_for_sg.get(sg, ""), cr_for_sg.get(sg, ""),
                str(a.get("top_umis", 0)), str(a.get("top_reads", 0)),
                str(a["total_umis"]), str(a["total_reads"]),
                f"{a.get('dominant_frac', 0):.4f}", str(a.get("n_sgrnas", 0)),
                str(a.get("second_sgid", "")), str(a.get("second_umis", 0)),
                f"{a.get('binom_p', 1):.3g}",
            ]) + "\n")

    print(f"Total cells with any sgRNA call: {n_total:,}")
    if h5_barcodes:
        print(f"Of {len(h5_barcodes):,} h5 cells, {n_in_h5:,} ({100*n_in_h5/len(h5_barcodes):.1f}%) have an sgRNA call")
        print(f"  ...of which {n_confident_in_h5:,} ({100*n_confident_in_h5/len(h5_barcodes):.1f}%) are 'confident' tier")

    # Summary by sgRNA / gene (for sanity check)
    sg_counts = Counter(a.get("top_sgid") for a in assigns.values() if a.get("tier") == "confident")
    gene_counts = Counter(gene_for_sg.get(s, "?") for s in sg_counts.elements())
    print("\n=== Confident assignments — top 15 sgRNAs ===")
    for s, c in sg_counts.most_common(15):
        print(f"  {s}\t{c}")
    print("\n=== Confident assignments — top 15 genes ===")
    for g, c in gene_counts.most_common(15):
        print(f"  {g}\t{c}")


if __name__ == "__main__":
    main()
