#!/bin/bash
# Run the task on line <n> of tasks.tsv (line 1 is the header).
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/abf_sdy
IFS=$'\t' read -r S C SG SE F < <(sed -n "${1}p" ${A}/tasks.tsv)
${RUN_ONE:-bash ${A}/run_one.sh} "$S" "$C" "$SG" "$SE" "$F"
