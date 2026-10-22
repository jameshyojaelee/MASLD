# Assembly on SYNTHETIC tables (debug_fake_outputs.R), tiny bootstrap; result written inside this debug leaf.
DSET=${H80}/debug/${DBG_SETUP_LEAF}
FAKE_CSG_MEAN=6 ${ENV}/bin/Rscript ${H80}/debug/debug_fake_outputs.R ${DSET} ${OUT} > ${OUT}/fake.log 2>&1
mkdir ${OUT}/assemble
${ENV}/bin/Rscript ${H80}/scripts/04_assemble.R --setup ${OUT}/fake_setup --arm_a ${OUT}/fake_armA --arm_b ${OUT}/fake_armB_pc1 --readout ${OUT}/fake_readout --out ${OUT}/assemble --n_boot 50 --result ${OUT}/RESULT_FAKE.md --skip_preamble TRUE > ${OUT}/assemble.log 2>&1
