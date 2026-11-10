"""SS S3f/S3g/S3h step 0. Allele counts from the per-library S2 VCFs (bcftools mpileup -R s1b_targets.tsv, FORMAT/AD),
for every included library (units.tsv include == TRUE). Stage is not read. Identity sites are skipped (used by S2 only).
Writes:
  geno_records.tsv.gz  raw records at the DIG1 job 62, PANX1 job 71 and TNFSF10 tag positions (job 62b input: ref,
                       alts, AD, GT if present, INDEL flag), one row per record.
  ase_counts.tsv.gz    SNP records at PANX1, TNFSF10, S3h and S3h-null positions: reads for the expected REF and ALT
                       allele (other alleles not counted), as in jobs 71b and 72b.
Missing VCF or position without a record = no reads.
Usage: python -I s3_allele_counts.py --input <S2 dir> --targets <out_s0d_targets dir> --out <new dir>
"""
import argparse
import gzip
import hashlib
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--input", required=True)
ap.add_argument("--targets", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
IN, TG, OUT = Path(a.input), Path(a.targets), Path(a.out)
OUT.mkdir(parents=True, exist_ok=False)
for line in (TG / "SHA256SUMS").read_text().splitlines():
    h, f = line.split(maxsplit=1)
    assert hashlib.sha256(Path(f).read_bytes()).hexdigest() == h, f"sha256 mismatch: {f}"
T = pd.read_csv(TG / "s1b_targets_annotated.tsv", sep="\t")
T = T[T.role != "identity"]
GENO_ROLES = {"dig1_62", "panx1_71", "tnfsf10"}
ASE_ROLES = {"panx1_71", "tnfsf10", "s3h", "s3h_null"}
geno_pos = set(zip(T[T.role.isin(GENO_ROLES)].chrom, T[T.role.isin(GENO_ROLES)].pos))
ase = T[T.role.isin(ASE_ROLES)].drop_duplicates(["chrom", "pos"])
ase_allele = {(c, p): (r, x) for c, p, r, x in zip(ase.chrom, ase.pos, ase.ref, ase.alt)}
U = pd.read_csv(IN / "units.tsv", sep="\t")
U = U[U.include.astype(str).str.upper() == "TRUE"]
n_missing = 0
with gzip.open(OUT / "geno_records.tsv.gz", "wt") as g, gzip.open(OUT / "ase_counts.tsv.gz", "wt") as s:
    g.write("sample_alias\tchrom\tpos\tref\talts\tad\tgt\tindel\n")
    s.write("sample_alias\tchrom\tpos\tnref\tnalt\n")
    for sa in U.sample_alias.astype(str):
        f = IN / "genotypes" / f"{sa}.vcf.gz"
        if not f.exists():
            n_missing += 1
            continue
        with gzip.open(f, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                key = (c[0], int(c[1]))
                if key not in geno_pos and key not in ase_allele:
                    continue
                fmt = dict(zip(c[8].split(":"), c[9].split(":")))
                ad = fmt.get("AD", "0")
                indel = "INDEL" in c[7].split(";")
                if key in geno_pos:
                    g.write(f"{sa}\t{c[0]}\t{c[1]}\t{c[3]}\t{c[4]}\t{ad}\t{fmt.get('GT', '.')}\t{int(indel)}\n")
                if key in ase_allele and not indel:
                    eref, ealt = ase_allele[key]
                    if c[3] != eref:
                        continue
                    v = [int(x) if x != "." else 0 for x in ad.split(",")]
                    al = c[4].split(",")
                    nalt = v[al.index(ealt) + 1] if ealt in al else 0
                    s.write(f"{sa}\t{c[0]}\t{c[1]}\t{v[0]}\t{nalt}\n")
(OUT / "summary.txt").write_text(f"included libraries {len(U)}; VCF missing {n_missing}\n")
print(f"included libraries {len(U)}; VCF missing {n_missing}")
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
