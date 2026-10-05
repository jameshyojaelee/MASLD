#!/usr/bin/env python3
"""Recheck small saved pilot assignments with current molecule restrictions."""
import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import groupby
import json
import os
from pathlib import Path

from count_verified_fragments import molecule_assignment


def main():
    import pysam
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pilot', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Requires allocated compute node')
    original = json.loads((args.pilot / 'DIAGNOSTIC_NOT_PRIMARY.json').read_text())
    report = {'diagnostic_only_not_primary_counts': True, 'source_BAM_rescanned': False,
              'sample_id': original['sample_id'], 'original_report': str(args.pilot / 'DIAGNOSTIC_NOT_PRIMARY.json'),
              'current_helper_sha256': hashlib.sha256(Path(__file__).with_name('count_verified_fragments.py').read_bytes()).hexdigest(),
              'comparisons': {}}
    for mode in ['baseline', 'multihit']:
        dispositions = Counter(); counts = Counter()
        with pysam.AlignmentFile(str(args.pilot / mode / 'assignments.namesorted.bam'), 'rb') as bam:
            for _, group in groupby(bam.fetch(until_eof=True), key=lambda r: r.query_name):
                groups = defaultdict(list)
                for r in group:
                    tags = dict(r.tags)
                    groups[tags.get('RG', '')].append({'NH': tags.get('NH'), 'HI': tags.get('HI'),
                        'XS': tags.get('XS'), 'XN': tags.get('XN'), 'XT': tags.get('XT'),
                        'unmapped': r.is_unmapped, 'supplementary': r.is_supplementary,
                        'chimeric': r.has_tag('SA') or r.has_tag('ch'), 'paired': r.is_paired,
                        'mate_unmapped': r.mate_is_unmapped, 'reference_id': r.reference_id,
                        'reference_start': r.reference_start, 'next_reference_id': r.next_reference_id,
                        'next_reference_start': r.next_reference_start,
                        'mate': 1 if r.is_read1 else (2 if r.is_read2 else 0)})
                for records in groups.values():
                    gene, reason = molecule_assignment(records, original['layout'] == 'PAIRED')
                    dispositions[reason] += 1
                    if gene:
                        counts[gene] += 1
        report['comparisons'][mode] = {'dispositions': dict(dispositions), 'accepted_all_sum': sum(counts.values())}
    with args.output.open('x') as handle:
        handle.write(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
