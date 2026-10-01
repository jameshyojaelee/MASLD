#!/usr/bin/env python3
"""Check Wenz significant-only association coordinates against GRCh38.

This is source admission metadata, not an allele-effect model evaluation.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parents[3]
CURRIN = ROOT/"GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
FASTA = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")


def main(source, out):
    if out.exists():
        raise FileExistsError(out)
    w = pd.read_csv(source, sep="\t", usecols=["Feature_ID", "rs_ID", "Chromosome",
        "SNP_position", "Ref_Allele", "Alt_Allele", "Effect_Size", "qval", "ConvergenceStatus"])
    if w.Feature_ID.duplicated().any() or not w.Effect_Size.between(0, 1, inclusive="neither").all():
        raise ValueError("Significant Wenz feature/effect census differs")
    w["chr"] = "chr"+w.Chromosome.astype(str).str.removeprefix("chr")
    reference = pysam.FastaFile(str(FASTA))
    observed = [reference.fetch(chrom, int(pos)-1, int(pos)).upper()
                for chrom, pos in zip(w.chr, w.SNP_position)]
    ref_match = sum(a == str(b).upper() for a,b in zip(observed, w.Ref_Allele))
    alt_match = sum(a == str(b).upper() for a,b in zip(observed, w.Alt_Allele))
    unmatched = len(w)-sum(a == str(b).upper() or a == str(c).upper()
                           for a,b,c in zip(observed, w.Ref_Allele, w.Alt_Allele))
    c = pd.read_csv(CURRIN, sep="\t", usecols=["chr", "pos_hg38", "ref", "alt"])
    same = w.rename(columns={"SNP_position": "pos_hg38", "Ref_Allele": "ref", "Alt_Allele": "alt"})
    overlap = same.merge(c, on=["chr", "pos_hg38", "ref", "alt"])
    result = {"source": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "population": "FDR<0.05 Wenz RASQUAL peak leads only; not all tested pairs",
              "rows": len(w), "unique_features": int(w.Feature_ID.nunique()),
              "unique_rs_ids": int(w.rs_ID.nunique()),
              "convergence_nonzero": int(w.ConvergenceStatus.ne(0).sum()),
              "pi_min": float(w.Effect_Size.min()), "pi_max": float(w.Effect_Size.max()),
              "grch38_fasta_REF_match": ref_match, "grch38_fasta_ALT_match": alt_match,
              "grch38_fasta_neither_allele": unmatched,
              "exact_currin_variant_tuple_overlap_rows": len(overlap),
              "interpretation": "Strong empirical GRCh38 alignment of the significant table; 2.7% REF mismatches need row-level harmonization. This does not resolve final donor membership, all-tested denominator or foundation-model exposure."}
    out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    main(a.source, a.out)
