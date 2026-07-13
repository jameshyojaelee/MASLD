#!/usr/bin/env python3
"""Convert STAR SJ.out.tab to regtools-compatible BED12 using unique reads."""
import argparse
p=argparse.ArgumentParser();p.add_argument("src");p.add_argument("dst");a=p.parse_args()
n=reads=0
with open(a.src) as src, open(a.dst+".tmp","w") as dst:
  for line in src:
    x=line.rstrip().split("\t")
    if len(x)<9: continue
    chrom,start,end,strand=x[0],int(x[1]),int(x[2]),int(x[3]); count=int(x[6])
    if count<=0 or strand not in (1,2): continue
    # LeafCutter's BED12 parser recovers intron (start-1,end+1). Eight-base
    # artificial anchors are sufficient because only intron coordinates/counts
    # enter clustering. STAR unique-read count is column 7.
    A=start-1; B=end+1; anchor=8; bed_start=A-anchor; bed_end=B+anchor-1
    if bed_start<0: continue
    st="+" if strand==1 else "-"; second=bed_end-bed_start-anchor
    dst.write(f"{chrom}\t{bed_start}\t{bed_end}\tSJ{n:08d}\t{count}\t{st}\t{bed_start}\t{bed_end}\t255,0,0\t2\t{anchor},{anchor}\t0,{second}\n")
    n+=1;reads+=count
import os;os.replace(a.dst+".tmp",a.dst)
print(f"junctions={n} unique_reads={reads}")
