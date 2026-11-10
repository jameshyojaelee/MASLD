# Trailing-zero trap check: restored lead position vs credible-set member span
import pandas as pd
A="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
c=pd.read_csv(A+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv",sep="\t",dtype={"locus":str})
d_lo=c.member_min_position-c.pos; d_hi=c.pos-c.member_max_position
far=(c.pos<c.member_min_position-1_000_000)|(c.pos>c.member_max_position+1_000_000)
print("sets",len(c),"restored",(c.trailing_zeros_restored>0).sum(),"pos >1 Mb outside member span",far.sum())
print("pos outside member span at all",((c.pos<c.member_min_position)|(c.pos>c.member_max_position)).sum())
# does locus string decode to pos?
dec=c.locus.str.split(".",n=1).str[1].astype(int)
print("locus decimal != restored pos",(dec!=c.pos).sum())
