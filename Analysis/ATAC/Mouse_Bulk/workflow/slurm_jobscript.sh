#!/bin/bash
# Snakemake SLURM jobscript — loads binary-only modules before each sub-job
#
# IMPORTANT: Python-based tools (MACS2, deepTools) are loaded via module
# commands INSIDE each rule's shell block to avoid PYTHONPATH conflicts
# with Snakemake's conda Python. Only binary tools loaded here.

# Load binary-only tools (no Python contamination)
module purge 2>/dev/null
module load Bowtie2/2.5.4-linux-x86_64
module load SAMtools/1.21
module load BEDTools/2.31.0-GCC-12.3.0
module load picard/3.0.0-Java-17
module load Subread/2.0.4-GCC-11.3.0

# Activate conda env for Snakemake post-processing
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Clear any PYTHONPATH to keep Snakemake's Python clean
unset PYTHONPATH

# Add fastp to PATH
export PATH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Mouse_Bulk/bin:$PATH"

# Execute the Snakemake-generated commands
{exec_job}
