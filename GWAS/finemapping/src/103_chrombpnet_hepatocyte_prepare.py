#!/usr/bin/env python3
"""Pool adult-hepatocyte fragments and convert peak BED3 to narrowPeak."""
import gzip, os
ROOT=os.environ.get("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATAC=os.path.join(ROOT,"Analysis/ATAC/Human_Multiome")
OUT=os.path.join(ROOT,"GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_fold0")
os.makedirs(OUT,exist_ok=True)
pooled=os.path.join(OUT,"hepatocyte.fragments.tsv.gz")
kept=total=0
with gzip.open(pooled+".tmp","wt") as dst:
  for donor in [f"D{i:02d}" for i in range(1,19)]:
    allow=set(open(os.path.join(ROOT,"GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/barcodes",donor+".hepatocyte.txt")).read().split())
    src=os.path.join(ATAC,"results/fragments",donor+"_fragments.tsv.gz")
    with gzip.open(src,"rt") as f:
      for line in f:
        total+=1; x=line.split("\t")
        if len(x)>=5 and x[3] in allow: dst.write(line); kept+=1
os.replace(pooled+".tmp",pooled)
if kept==0: raise SystemExit("no fragments matched hepatocyte barcode allowlists")
bed3=os.path.join(ATAC,"results/label_transfer/cell_type_peak_sets_v2/Hepatocytes_peaks.bed")
np=os.path.join(OUT,"hepatocyte.peaks.narrowPeak")
with open(bed3) as src, open(np,"w") as dst:
  for i,line in enumerate(src,1):
    c,s,e=line.rstrip().split("\t")[:3]; s=int(s);e=int(e); summit=(e-s)//2
    dst.write(f"{c}\t{s}\t{e}\thep_peak_{i}\t1000\t.\t10\t-1\t-1\t{summit}\n")
fai="/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa.fai"
with open(fai) as src, open(os.path.join(OUT,"GRCh38.chrom.sizes"),"w") as dst:
  for line in src:
    x=line.split("\t");
    # Coverage construction sees alternate contigs present in fragments. Keep the
    # complete FASTA index here; the fold JSON still limits model fitting/evaluation
    # to explicitly assigned primary chromosomes.
    dst.write(x[0]+"\t"+x[1]+"\n")
with open(os.path.join(OUT,"prep_qc.tsv"),"w") as f:f.write(f"total_fragments\t{total}\nhepatocyte_fragments\t{kept}\nretained_fraction\t{kept/total:.6f}\n")
print(f"kept {kept}/{total} fragments")
