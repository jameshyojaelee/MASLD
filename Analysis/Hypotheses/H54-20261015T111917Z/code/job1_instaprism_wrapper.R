# H54 job 1, step 2. Runs the unchanged InstaPrism script with set.seed(20260923)
# immediately before it (prespec "Seed"). Called as:
#   bin/R --no-save --no-restore --slave --file=job1_instaprism_wrapper.R --args <ref> <counts> <outdir> <dataset> celltype
# commandArgs(trailingOnly = TRUE) inside the sourced script sees the same --args.
set.seed(20260923L)
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Deconvolution/scripts/12_run_bayesprism.R")
