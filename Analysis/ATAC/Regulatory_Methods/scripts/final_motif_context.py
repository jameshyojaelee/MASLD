#!/usr/bin/env python3
"""Reuse exact-peak frozen motif scans for descriptive all-donor link maps."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path

import annotate_motif_families as families
from final_annotation_paths import annotation_paths

A=Path(__file__).resolve().parents[1]
ROOT=A.parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--cohort',required=True)
    p.add_argument('--lineage',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    args.output.mkdir(parents=True,exist_ok=False)
    directory=annotation_paths(args.candidate)['disease_links']/f'{args.cohort}_{args.lineage}'
    assert json.loads((directory/'completion.json').read_text())['status']=='completed'
    source=args.output/'source'; source.mkdir()
    (source/'completion.json').write_text(json.dumps(dict(status='completed',pilot=False,
        interpretation='all_donor_descriptive_accessibility_links_not_prediction_folds')))
    target=source/f'{args.cohort}_{args.lineage}';target.mkdir()
    columns=['outer_fold','gene_id','peak_id','method','rank','score','score_unit','evidence_state','cohort','lineage']
    with gzip.open(directory/'links.tsv.gz','rt') as inp,(target/'links.tsv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=columns,delimiter='\t'); writer.writeheader()
        for row in csv.DictReader(inp,delimiter='\t'):
            row.update(outer_fold='all_donors_descriptive',evidence_state='candidate_accessibility_link')
            writer.writerow({key:row[key] for key in columns})
    motif=ROOT/'Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/chromvar'
    families.annotate([source],motif,args.output/'annotation')
    (args.output/'completion.json').write_text(json.dumps(dict(status='completed',
        interpretation='TF_family_compatibility_only_no_TF_identity_activity_binding_or_causality',scientific_adoption=False),indent=2))


if __name__=='__main__': main()
