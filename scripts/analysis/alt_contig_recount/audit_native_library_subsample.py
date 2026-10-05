#!/usr/bin/env python3
"""Diagnostic native-GTF molecule audit; never produces primary counts.

Scan the entire source BAM and retain every record of a deterministic 1/1000
QNAME/RG sample. Region restriction and alignment-level decimation are forbidden.
"""
import argparse
from collections import Counter, defaultdict
import csv
from functools import lru_cache
import hashlib
from itertools import groupby
import json
import os
from pathlib import Path
import subprocess

from count_verified_fragments import molecule_assignment

BIN = Path('/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin')
SEED = 'native-library-diagnostic-20261005'


@lru_cache(maxsize=250000)
def selected(qname, read_group, seed=SEED):
    key = json.dumps([seed, qname, read_group], ensure_ascii=False, separators=(',', ':'))
    return int.from_bytes(hashlib.sha256(key.encode()).digest(), 'big') % 1000 == 0


def run(command, log):
    with log.open('a') as handle:
        handle.write(json.dumps([str(x) for x in command]) + '\n')
        handle.flush()
        subprocess.run([str(x) for x in command], stdout=handle, stderr=handle, check=True)


def stat(path):
    p = Path(path).resolve(); s = p.stat()
    return {'path': str(p), 'size_bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}


def main():
    import pysam
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--native-gtf', type=Path, required=True)
    parser.add_argument('--expected-gtf-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dispatch-worker', action='store_true',
                        help='Execute a pilot explicitly rerouted from the serial queue')
    parser.add_argument('--independent-worker', action='store_true',
                        help='Execute a pilot rerouted from the external serial queue')
    args = parser.parse_args()
    independent_registry = args.output.parent / 'independent_dispatch.json'
    if independent_registry.exists() and args.sample_id in json.loads(independent_registry.read_text())['sample_ids'] and not args.independent_worker:
        print('SKIP_QUEUE_INDEPENDENTLY_DISPATCHED', args.sample_id, flush=True)
        return
    registry = args.output.parent / 'external_dispatch.json'
    if registry.exists() and args.sample_id in json.loads(registry.read_text())['sample_ids'] and not args.dispatch_worker:
        print('SKIP_SERIAL_EXTERNALLY_DISPATCHED', args.sample_id, flush=True)
        return
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Requires existing allocated compute node')
    with args.manifest.open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter='\t') if r['sample_id'] == args.sample_id]
    if len(rows) != 1:
        raise ValueError('Manifest must uniquely identify sample')
    row = rows[0]
    sampler_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    helper_sha256 = hashlib.sha256(Path(__file__).with_name('count_verified_fragments.py').read_bytes()).hexdigest()
    source = Path(row['bam_path'])
    layout = row['layout'].upper()
    if layout not in {'SINGLE', 'PAIRED'}:
        raise ValueError('Expected SINGLE or PAIRED layout')
    strand = int(row['strandedness'])
    if strand not in {0, 1, 2}:
        raise ValueError('Invalid native featureCounts strand rule')
    gtf_hash = hashlib.sha256(args.native_gtf.read_bytes()).hexdigest()
    if gtf_hash != args.expected_gtf_sha256:
        raise ValueError('Native GTF differs from frozen annotation resolver input')
    args.output.mkdir(parents=True, exist_ok=False)
    sample_bam = args.output / 'whole_source_hash_sample.bam'
    stats = Counter()
    before = stat(source)
    with pysam.AlignmentFile(str(source), 'rb') as bam:
        with pysam.AlignmentFile(str(sample_bam), 'wb', template=bam) as out:
            for record in bam.fetch(until_eof=True):
                stats['source_records_scanned'] += 1
                if stats['source_records_scanned'] % 1000000 == 0:
                    print(args.sample_id, 'scanned_records', stats['source_records_scanned'],
                          'retained_records', stats['sample_records'], flush=True)
                rg = record.get_tag('RG') if record.has_tag('RG') else ''
                if selected(record.query_name, rg):
                    out.write(record)
                    stats['sample_records'] += 1
                    for label, flag in [('unmapped', record.is_unmapped), ('secondary', record.is_secondary),
                                        ('supplementary', record.is_supplementary), ('SA_tag', record.has_tag('SA')),
                                        ('ch_tag', record.has_tag('ch')), ('missing_NH', not record.has_tag('NH')),
                                        ('missing_HI', not record.has_tag('HI')),
                                        ('mate_unmapped', record.is_paired and record.mate_is_unmapped),
                                        ('cross_contig_mate', record.is_paired and not record.mate_is_unmapped
                                         and record.reference_id != record.next_reference_id),
                                        ('improper_pair', record.is_paired and not record.is_proper_pair)]:
                        stats[label] += int(flag)
    if before != stat(source):
        raise RuntimeError('Source BAM changed during sampling')
    log = args.output / 'commands.log'
    flags = ['-T', '1', '-s', str(strand), '-a', args.native_gtf]
    if layout == 'PAIRED':
        flags += ['-p', '--countReadPairs', '-B']
    comparisons = {}
    for mode in ['baseline', 'multihit']:
        folder = args.output / mode; folder.mkdir()
        # Separate BAM names prevent featureCounts assignment output collisions.
        bam_link = folder / 'sample.bam'; bam_link.symlink_to(sample_bam.resolve())
        count_path = folder / 'featurecounts.tsv'
        command = [BIN / 'featureCounts', *flags, '-R', 'BAM', '-o', count_path]
        if mode == 'multihit':
            command += ['-M']
        run([*command, bam_link], log)
        assigned = folder / 'sample.bam.featureCounts.bam'
        sorted_bam = folder / 'assignments.namesorted.bam'
        run([BIN / 'samtools', 'sort', '-N', '-@', '1', '-m', '512M', '-o', sorted_bam, assigned], log)
        dispositions = Counter(); counts = Counter(); unique = Counter()
        with pysam.AlignmentFile(str(sorted_bam), 'rb') as bam:
            for _, group in groupby(bam.fetch(until_eof=True), key=lambda r: r.query_name):
                groups = defaultdict(list)
                for record in group:
                    tags = dict(record.tags)
                    groups[tags.get('RG', '')].append({'NH': tags.get('NH'), 'HI': tags.get('HI'),
                        'XS': tags.get('XS'), 'XN': tags.get('XN'), 'XT': tags.get('XT'),
                        'unmapped': record.is_unmapped, 'supplementary': record.is_supplementary,
                        'chimeric': record.has_tag('SA') or record.has_tag('ch'), 'paired': record.is_paired,
                        'mate_unmapped': record.mate_is_unmapped,
                        'reference_id': record.reference_id, 'reference_start': record.reference_start,
                        'next_reference_id': record.next_reference_id, 'next_reference_start': record.next_reference_start,
                        'mate': 1 if record.is_read1 else (2 if record.is_read2 else 0)})
                for records in groups.values():
                    gene, reason = molecule_assignment(records, layout == 'PAIRED')
                    dispositions[reason] += 1
                    if gene:
                        counts[gene] += 1
                        if reason == 'accepted_unique':
                            unique[gene] += 1
        fc = {}
        with count_path.open() as handle:
            for r in csv.DictReader((line for line in handle if not line.startswith('#')), delimiter='\t'):
                fc[r['Geneid']] = int(r[next(reversed(r))])
        with (folder / 'diagnostic_gene_comparison.tsv').open('x') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['gene_id', 'native_featureCounts', 'accepted_unique_molecules',
                             'all_accepted_molecules', 'unique_minus_featureCounts'])
            for gene in sorted(set(fc) | set(counts)):
                writer.writerow([gene, fc.get(gene, 0), unique[gene], counts[gene], unique[gene] - fc.get(gene, 0)])
        comparisons[mode] = {'dispositions': dict(dispositions), 'native_featureCounts_sum': sum(fc.values()),
                             'accepted_unique_sum': sum(unique.values()), 'accepted_all_sum': sum(counts.values()),
                             'gene_rows_with_unique_count_difference': sum(unique[g] != fc.get(g, 0) for g in set(fc) | set(unique))}
    versions = {}
    for binary, flag in [('samtools', '--version'), ('featureCounts', '-v')]:
        result = subprocess.run([str(BIN / binary), flag], capture_output=True, text=True, check=True)
        versions[binary] = (result.stdout + result.stderr).strip()
    report = {'diagnostic_only_not_primary_counts': True, 'sample_id': args.sample_id,
              'seed': SEED, 'hash_modulus': 1000, 'sampling_unit': 'QNAME_plus_RG_all_source_records',
              'source': before, 'manifest': stat(args.manifest), 'native_gtf': stat(args.native_gtf),
              'native_gtf_sha256': gtf_hash,
              'molecule_assignment_source_sha256_at_start': helper_sha256,
              'sampler_source_sha256_at_start': sampler_sha256,
              'sampler_selection_cache_maxsize': 250000,
              'layout': layout, 'strandedness': strand, 'versions': versions, 'pysam': pysam.__version__,
              'stats': dict(stats), 'comparisons': comparisons,
              'difference_interpretation': 'Native unresolved GTF and complete-alignment acceptance differ; dispositions describe rejection reasons, without resolving ALT equivalence or certifying counts.'}
    (args.output / 'DIAGNOSTIC_NOT_PRIMARY.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
