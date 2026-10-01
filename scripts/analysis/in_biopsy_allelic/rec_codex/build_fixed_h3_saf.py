"""Export the frozen complete H3 region axis as unstranded-count SAF metadata.

No molecular counts, alignment files, receiving values or prediction fitting.
The dummy positive strand is ignored only with featureCounts -s 0.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
AXIS = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/h3k27ac_feature_axis.tsv'
GENOME = REC / 'ensembl98_primary_genome_21998416'
FAI = GENOME / 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
GUARDS = {
    AXIS: 'cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986',
    FAI: '57f6ed6f8b07437708fff09814489cfda4b976652afabe2c9b1d54463bcba0ad',
    GENOME / 'summary.json': '9f93a8291cf676bbf1bd1bc552d2097b3e0f9850a05f86431e1d93df65392f6b',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID'), 'Compute allocation required')
    require(not args.output.exists(), 'Refusing output overwrite')
    hashes = {str(p): sha(p) for p in GUARDS}
    require(all(hashes[str(p)] == expected for p, expected in GUARDS.items()), 'Frozen metadata changed')
    genome = json.loads((GENOME / 'summary.json').read_text())
    # The genome sequence was already hashed and its FAI generated on compute.
    require(genome['fasta_sha256'] == '78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c'
            and genome['faidx_sha256'] == GUARDS[FAI], 'Wrong declared primary genome')
    contigs = {}
    for line in FAI.read_text().splitlines():
        name, length, *_ = line.split('\t')
        require(name not in contigs and int(length) > 0, 'Invalid genome index')
        contigs[name] = int(length)
    require(len(contigs) == 194, 'Genome contig census changed')
    require(contigs == genome['contig_lengths'], 'Genome summary and FAI differ')
    with AXIS.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    require(len(rows) == 96460, 'Complete target region census changed')
    keys, saf, crosswalk, chromosome_counts = set(), [], [], {}
    for i, row in enumerate(rows):
        key = row['opaque_source_feature_key']
        require(int(row['h3k27ac_feature_index']) == i and key not in keys, 'Region order or uniqueness changed')
        keys.add(key)
        match = re.fullmatch(r'chr([0-9]+|X|Y):([0-9]+)-([0-9]+)', key)
        require(match is not None, 'Unsupported original region key: ' + key)
        chromosome, start, end = match.groups()
        start, end = int(start), int(end)
        require(chromosome in contigs and 1 <= start <= end <= contigs[chromosome], 'Region outside declared genome')
        require(chromosome in [str(c) for c in range(1, 23)] + ['X', 'Y'], 'Unexpected chromosome')
        # SAF and original keys are both one-based inclusive. Do not shift.
        saf.append(dict(GeneID=key, Chr=chromosome, Start=start, End=end, Strand='+'))
        crosswalk.append(dict(h3k27ac_feature_index=i, opaque_source_feature_key=key,
                              alignment_contig=chromosome, start1=start, end1=end,
                              bed_start0=start - 1, bed_end0=end, interval_width_bp=end-start+1))
        chromosome_counts[chromosome] = chromosome_counts.get(chromosome, 0) + 1
    require(sum(chromosome_counts.get(str(c), 0) for c in range(1, 23)) == 94026
            and chromosome_counts['X'] == 2352 and chromosome_counts['Y'] == 82,
            'Frozen chromosome census changed')
    args.output.mkdir(parents=True, exist_ok=False)
    for name, records in [('fixed_h3_regions.saf', saf), ('region_coordinate_crosswalk.tsv', crosswalk)]:
        with (args.output / name).open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(records)
    (args.output / 'executed_source.py').write_bytes(Path(__file__).read_bytes())
    summary = dict(schema_version='complete-fixed-h3-saf-v1', status='fixed_annotation_complete',
                   regions=len(saf), contigs=len(contigs), chromosome_region_counts=chromosome_counts,
                   coordinates='one-based inclusive original region keys and SAF; BED conversion start0=start1-1,end0=end1',
                   strand='dummy +; biological strand not inferred; count only with -s 0',
                   prospective_count_unit='single primary nonsupplementary mapped read overlapping one fixed region',
                   prospective_count_qualifications='MAPQ>=10; no -O or -M; no paired-fragment counting, receiving peak calls, Input subtraction or PCR deduplication',
                   molecular_values_read=False, alignments_read=False, receiving_values_read=False,
                   prediction_model_retrained=False, target_definition_changed=False,
                   reference_measurement_claim='annotation adapter only; not proof of original CUT&RUN-count reconstruction or external ChIP validity',
                   input_sha256=hashes, script_sha256=sha(Path(__file__)), python=sys.version,
                   slurm_job_id=os.environ['SLURM_JOB_ID'],
                   output_sha256={p.name: sha(p) for p in args.output.iterdir() if p.is_file()})
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'status': summary['status'], 'regions': len(saf)}))


if __name__ == '__main__':
    main()
