#!/bin/bash
# Orchestrator for the CIRRHOSIS-EXCLUDED hepatocyte Hotspot re-run (remediation R1).
# Chains:  501 (--exclude-stage Cirrhosis, bigmem/interactive/512G)
#       -> 505b phenotype correlation (cpu/32G)  [afterok]
#       -> 506b bulk replication vs C2 canonical  (same cpu job, sequential)
#
# 501 is heavy (657k hepatocytes need ~512G > cpu's 210G cap), so it uses
# plain bigmem (nslab QOS) mirroring the canonical 501_run_hotspot.sbatch. The
# light R steps run together in one cpu job, gated on 501 via afterok.
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
DIR=Analysis/SingleCell/scripts/hotspot_modules

JID_501=$(sbatch --parsable "$DIR/501_nocirrhosis_hep.sbatch")
echo "Submitted 501 (nocirrhosis hepatocytes, bigmem/nslab/512G): $JID_501"

JID_DOWN=$(sbatch --parsable --dependency=afterok:"$JID_501" "$DIR/505b_506b_nocirr.sbatch")
echo "Submitted 505b+506b (cpu, afterok:$JID_501): $JID_DOWN"

echo ""
echo "Dependency chain:  $JID_501 (501)  ->  $JID_DOWN (505b -> 506b)"
echo "Track with:  squeue --name=hotspot -u \$USER"
