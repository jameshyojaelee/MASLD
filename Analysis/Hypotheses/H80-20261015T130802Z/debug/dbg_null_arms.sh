# Null-only debug of jobs 1-3 against the debug setup (lambda = 0 everywhere; no outcome computed).
DSET=${H80}/debug/${DBG_SETUP_LEAF}
ROWS=$(${ENV}/bin/Rscript -e 'suppressMessages(library(data.table)); w <- fread(commandArgs(TRUE)[1]); o <- w[order(nE)]; cat(paste(unique(c(o$row[1], o$row[ceiling(nrow(o)/2)], o$row[nrow(o)])), collapse=","))' ${DSET}/windows.tsv)
BROW=$(${ENV}/bin/Rscript -e 'suppressMessages(library(data.table)); w <- fread(commandArgs(TRUE)[1]); cat(w[armB == TRUE][order(nE)][ceiling(.N/2), row])' ${DSET}/windows.tsv)
echo "Arm A debug rows ${ROWS}; Arm B row ${BROW}"
mkdir ${OUT}/armA ${OUT}/armB ${OUT}/readout
/usr/bin/time -v ${ENV}/bin/Rscript ${H80}/scripts/01_arm_a.R --setup ${DSET} --out ${OUT}/armA --workers 3 --only_rows ${ROWS} --n_rep 2 --null_only TRUE --skip_preamble TRUE > ${OUT}/armA.log 2>&1
/usr/bin/time -v ${ENV}/bin/Rscript ${H80}/scripts/02_arm_b_pc1.R --setup ${DSET} --out ${OUT}/armB --workers 2 --only_b ${BROW} --only_p 1 --n_rep 2 --null_only TRUE --skip_preamble TRUE > ${OUT}/armB.log 2>&1
/usr/bin/time -v ${ENV}/bin/Rscript ${H80}/scripts/03_lead_z.R --setup ${DSET} --out ${OUT}/readout --workers 4 --only_chr 22 --skip_preamble TRUE > ${OUT}/readout.log 2>&1
