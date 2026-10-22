#!/bin/bash
# H74 prespec: group counts from the NG13 gate per-locus table (no outcome). Login-node awk on an 8.7 kB file.
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG13
awk -F'\t' 'NR>1{n++; if($7==1)a++; else b++; if($8==0)z++; if($5==$7)u++; if($6<$5)p++; if($4>$3)m++} END{print "loci",n,"n_h4==1",a,"n_h4>=2",b,"n_h3==0",z,"all_tested_coloc",u,"loci_with_nonPC_tested",p,"multi_snp_clusters",m}' locus_structure.tsv
