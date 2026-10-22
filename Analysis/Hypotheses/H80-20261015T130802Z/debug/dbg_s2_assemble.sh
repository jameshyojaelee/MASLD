# Assembly on the REAL setup leaf (STATUS STOP:S2) with skipped Arm A / Arm B leaves, to confirm job 4 writes the
# STOPPED verdict without error. Nothing simulated exists; result written inside this debug leaf only.
mkdir ${OUT}/armA_skipped ${OUT}/armB_skipped ${OUT}/assemble
echo "SKIPPED: STOP:S2:PC-2_eqtl_fit_reproduction_failed" > ${OUT}/armA_skipped/SKIPPED
echo "SKIPPED: STOP:S2:PC-2_eqtl_fit_reproduction_failed" > ${OUT}/armB_skipped/SKIPPED
${ENV}/bin/Rscript ${H80}/scripts/04_assemble.R --setup ${H80}/run_setup --arm_a ${OUT}/armA_skipped --arm_b ${OUT}/armB_skipped --readout ${H80}/run_readout --out ${OUT}/assemble --n_boot 50 --result ${OUT}/RESULT_S2TEST.md --skip_preamble TRUE > ${OUT}/assemble.log 2>&1
