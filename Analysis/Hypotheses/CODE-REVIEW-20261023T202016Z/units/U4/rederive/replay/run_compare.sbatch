#!/bin/bash
#SBATCH --job-name=u4_replay_cmp
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/logs/%j.out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/replay
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla $R/compare_replay.R
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla -e 'library(data.table); options(width=220); x<-fread("'$R'/replay_comparison.tsv"); print(x[,.(n=.N, ok=sum(status=="ok"), calls_same=sum(call_same,na.rm=TRUE), method_same=sum(method_same,na.rm=TRUE), na_same=sum(na_pattern_same,na.rm=TRUE), max_pp_diff=suppressWarnings(max(max_abs_pp_diff[is.finite(max_abs_pp_diff)])), top_snp_same=sum(top_snp_release==top_snp_replay,na.rm=TRUE), ncs_same=sum(n_cs_pairs_release==n_cs_pairs_replay | (is.na(n_cs_pairs_release)&is.na(n_cs_pairs_replay))), elig_same=sum(n_cs_pairs_eligible_release==n_cs_pairs_eligible_replay | (is.na(n_cs_pairs_eligible_release)&is.na(n_cs_pairs_eligible_replay)),na.rm=TRUE), nsnp_same=sum(n_snps_release==n_snps_replay,na.rm=TRUE)), by=stratum]); print(x[call_same==FALSE | method_same==FALSE | na_pattern_same==FALSE | max_abs_pp_diff>1e-6 | status!="ok", .(task,stratum,gene,status,max_abs_pp_diff,method_release,method_replay)])'
