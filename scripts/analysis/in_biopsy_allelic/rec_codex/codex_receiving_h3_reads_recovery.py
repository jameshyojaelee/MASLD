"""Recover fixed H3 processing; preserve original measurements and failed files.

Reused count payloads are hashed as bytes on compute, never parsed or copied.
No predictions, receiving count values, or accuracy are deserialized here.
"""
import argparse
from http.client import IncompleteRead
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import traceback
from urllib.error import HTTPError, URLError

import codex_receiving_h3_reads as helper

HERE = Path(__file__).resolve().parent
REC = helper.REC
ORIGINAL = REC / 'codex_receiving_h3_reads_21998882'
SOURCE_SHA = 'f641ea6485515b7f986f3ce279f449c5e6b85e99f676c3a5ac959fc93bf8b509'
ORIGINAL_GUARDS = {
    ORIGINAL/'protocol.json': 'd6394b1f4a938a1f5681b727287d4d7676e67af2f53c9329585edc62af19a609',
    ORIGINAL/'failure.json': '7f52a8cc08a4b272d4b79d0074d8055527c9d043cc54edcd026e0b6e7094d979',
    ORIGINAL/'codex_receiving_h3_reads.py': SOURCE_SHA,
}
REUSED = {
    'B1': '456fc548c536ef5b7d124287eed9b5eedb3b2061ce14643bd10204cfe2ddcc27',
    'B7': 'bc8d1893bf82c7ad044f0eeb380a187324cbfd6a075221d44b8c3aefd89aad2d',
    'B8': 'cdde39633d9d219cd6ded203cceedc2e56e763f083137c46975a51f26ced5998',
    'B9': '31f34bb915b5f78eac6c2211dac40ec7dc990cfba0b83fe58799671949b17e5a',
    'B10': 'cb670d9f32506c7598ed82f55dcef8c5286c6f1b02503b6540dcbf388785e35e',
    'B12': 'cd1ff9c3bfdb8a1e6fd980b8643ebc9a3abac7e9c44c9b85e47d321284d14c1b',
    'B15': 'be0ce8ad4ae700c7a2701c898b24471d9d5c0fc1d906756735df9edb359c9a44',
    'B19': '1b609114c285f5ea1953d9973b72e5739ebfb01226a113aa8c70387052f70467',
    'B21': '8c9e368c84547855733f039e47ef731f3260e287caae9e2eead12e8a714ca0fb',
    'B22': '0b527b6c622851fff4f98eb9c45422827a5100292b852d8a78a3123447bc8c61',
    'B24': '56d10c31f1a41161b0927afda16d64800aa4e0f70e53f92c10675c4f656610ee',
    'B26': '92249f420a9dea86bb099aa204bbaf0d2d3f7dfb33fdf48cdf501b84361326d0',
}
require = helper.require
sha = helper.sha
read_json = helper.read_json
guard = helper.guard
write_json = helper.write_json


def commands(row, base, bowtie, index):
    run_id = row['run_accession']
    fastq = base/'raw'/Path(row['fastq_url']).name
    aligned = base/(run_id+'.aligned.bam')
    filtered = base/(run_id+'.filtered.bam')
    align = [bowtie, '--end-to-end', '--sensitive', '--phred33', '--seed', str(helper.SEED),
             '-p', '7', '-x', index['index_prefix'], '-U', str(fastq),
             '--rg-id', run_id, '--rg', 'SM:'+row['donor']]
    to_bam = [str(helper.ENV/'samtools'), 'view', '-b', '-o', str(aligned), '-']
    prefilter = [str(helper.ENV/'samtools'), 'view', '-b', '-F', '2308', '-q', '10',
                 '-o', str(filtered), str(aligned)]
    return fastq, aligned, filtered, align, to_bam, prefilter


def count_command(base, bindings):
    return [str(helper.ENV/'featureCounts'), '-T', '8', '-F', 'SAF', '-a', str(helper.SAF),
            '-s', '0', '-Q', '10', '--primary', '-o', str(base/'native_featurecounts.tsv'), *bindings]


def admit_reused(frozen, rows, versions, bowtie, index):
    guard(ORIGINAL_GUARDS)
    protocol = read_json(ORIGINAL/'protocol.json')
    failure = read_json(ORIGINAL/'failure.json')
    require(not (ORIGINAL/'summary.json').exists() and failure['donor']=='B36'
            and failure['phase']=='acquisition' and tuple(failure['completed'])==tuple(REUSED)
            and failure['partials_preserved'], 'Original processing failure boundary changed')
    require(tuple(protocol['fixed_donors'])==helper.DONORS and protocol['h3_only_files']==202
            and protocol['compressed_bytes']==56783483922 and protocol['versions']==versions
            and protocol['alignment_seed']==helper.SEED and protocol['RNA_predictions_frozen_before_H3']
            and not protocol['model_fit'] and not protocol['predictive_accuracy_evaluated'],
            'Original processing cohort, versions, seed, or measurement boundary changed')
    bound = {Path(p):h for p,h in protocol['input_sha256'].items()}
    require(all(bound.get(p)==h for p,h in frozen.items()), 'Original prerequisite/binary binding changed')
    require(bound.get(HERE/'codex_receiving_h3_reads.py')==SOURCE_SHA,
            'Original executed processing source differs')
    guard(bound)
    reuse_guards = dict(ORIGINAL_GUARDS)
    receipts = {}
    for donor, expected in REUSED.items():
        base = ORIGINAL/donor
        receipt = base/'summary.json'
        require(sha(receipt)==expected, 'Successful original donor receipt changed')
        s = read_json(receipt)
        require(s['donor']==donor and s['h3_reads_processed'] and not s['predictions_refitted']
                and not s['count_axis_assembly_or_target_evaluation_run'], 'Original donor not complete')
        donor_rows = sorted((r for r in rows if r['donor']==donor), key=lambda r:r['run_accession'])
        bindings, filters, aligns = {}, {}, []
        require(len(s['native_read_files'])==len(donor_rows), 'Original acquisition coverage changed')
        for row, acquired in zip(donor_rows, s['native_read_files']):
            fastq, _, filtered, align, to_bam, prefilter = commands(row, base, bowtie, index)
            require(all(acquired.get(k)==v for k,v in row.items()) and acquired['path']==str(fastq)
                    and acquired['gzip_crc_verified'], 'Original acquisition identity/check boundary changed')
            bindings[str(filtered)] = {k:row[k] for k in
                ('run_accession', 'donor', 'GSM', 'experiment_accession', 'source_biosample')}
            filters[str(filtered)] = prefilter
            aligns.append(dict(alignment=align, to_bam=to_bam, prefilter=prefilter))
        require(s['columns']==bindings and s['prefilter_commands']==filters
                and s['alignment_commands']==aligns and s['count_command']==count_command(base, bindings)
                and s['count_table']==str(base/'native_featurecounts.tsv'), 'Original counting recipe changed')
        reuse_guards[receipt] = expected
        reuse_guards[base/'native_featurecounts.tsv'] = s['count_sha256']
        reuse_guards[base/'native_featurecounts.tsv.summary'] = s['count_summary_sha256']
        receipts[donor] = receipt
    # Byte verification only; no receiving value or BAM parsing.
    guard(reuse_guards)
    return reuse_guards, receipts


def acquire_with_retries(row, fastq, output, base, donor_cap):
    for attempt in range(1, 4):
        helper.bound(output, base, donor_cap)
        try:
            acquired = helper.acquire(row, fastq, output, base, donor_cap)
            acquired['acquisition_attempts'] = attempt
            return acquired
        except (ConnectionError, TimeoutError, IncompleteRead, URLError) as exc:
            # Refusal on HTTP status, identity, checksum, CRC, or storage errors is immediate.
            if isinstance(exc, HTTPError):
                raise
            partial = Path(str(fastq)+'.partial')
            retained = Path(str(fastq)+'.partial.attempt'+str(attempt))
            if partial.exists():
                require(not retained.exists(), 'Refuse retry partial overwrite')
                partial.rename(retained)
            write_json(base/(row['run_accession']+'_retry_'+str(attempt)+'.json'),
                       dict(run_accession=row['run_accession'], attempt=attempt,
                            error_type=type(exc).__name__, error=str(exc),
                            retained_partial=str(retained) if retained.exists() else None))
            helper.bound(output, base, donor_cap)
            if attempt==3:
                raise
            time.sleep(10*attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source-sha256', required=True)
    args = parser.parse_args()
    require(sha(Path(__file__))==args.expected_source_sha256, 'Recovery source differs from submitted source')
    require(os.environ.get('SLURM_JOB_ID', '').isdigit()
            and int(os.environ.get('SLURM_CPUS_PER_TASK', '0'))>=8, 'Require eight allocated CPUs')
    require(sha(HERE/'codex_receiving_h3_reads.py')==SOURCE_SHA, 'Frozen original helper changed')
    frozen, rows, index = helper.prerequisites()
    bowtie = shutil.which('bowtie2')
    require(bowtie, 'Load Bowtie2/2.5.4-linux-x86_64')
    versions = {}
    for name, cmd in [('bowtie2', [bowtie, '--version']),
                      ('samtools', [str(helper.ENV/'samtools'), '--version']),
                      ('featureCounts', [str(helper.ENV/'featureCounts'), '-v'])]:
        p = subprocess.run(cmd, capture_output=True, text=True, check=True)
        versions[name] = p.stdout+p.stderr
    frozen.update({Path(bowtie):sha(bowtie), helper.ENV/'samtools':sha(helper.ENV/'samtools'),
                   helper.ENV/'featureCounts':sha(helper.ENV/'featureCounts')})
    reuse_guards, receipts = admit_reused(frozen, rows, versions, bowtie, index)
    frozen.update({HERE/'codex_receiving_h3_reads.py':SOURCE_SHA,
                   HERE/'run_codex_receiving_h3_reads.sbatch':sha(HERE/'run_codex_receiving_h3_reads.sbatch'),
                   Path(__file__):args.expected_source_sha256,
                   HERE/'run_codex_receiving_h3_reads_recovery.sbatch':sha(HERE/'run_codex_receiving_h3_reads_recovery.sbatch')})
    output = REC/('codex_receiving_h3_reads_recovery_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    native_output_bytes = helper.output_bytes

    def combined_output_bytes(path):
        size = native_output_bytes(path)
        if path==output:
            size += native_output_bytes(ORIGINAL)
        elif path.parent==output and path.name in helper.DONORS and (ORIGINAL/path.name).exists():
            size += native_output_bytes(ORIGINAL/path.name)
        return size

    # Existing periodic caps now count original retained bytes plus new bytes.
    helper.output_bytes = combined_output_bytes
    donor_caps = {d:5*sum(int(r['fastq_bytes']) for r in rows if r['donor']==d)+2_000_000_000
                  for d in helper.DONORS}
    original_protocol = read_json(ORIGINAL/'protocol.json')
    protocol = dict(original_protocol)
    protocol.update(input_sha256={str(p):h for p,h in frozen.items()},
                    reused_input_sha256={str(p):h for p,h in reuse_guards.items()},
                    original_processing_directory=str(ORIGINAL), reused_donors=list(REUSED),
                    reacquired_donors=[d for d in helper.DONORS if d not in REUSED],
                    recovery_source_sha256=args.expected_source_sha256,
                    download_max_attempts=3, retry_delay_seconds=[10, 20],
                    retry_partial_policy='preserve each failed attempt; restart same complete URL without Range',
                    count_payloads_reused_in_place=True, count_payloads_parsed_or_copied=False,
                    storage_scope='all original retained files plus all new files; per-donor also combines old/new')
    write_json(output/'protocol.json', protocol)
    (output/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    write_json(output/'reused_receipts.json', dict(donor_receipt_paths={d:str(p) for d,p in receipts.items()},
               donor_summary_sha256=REUSED, original_source_sha256=SOURCE_SHA))
    completed = list(REUSED)
    donor = None
    phase = 'initialization'
    try:
        for donor in helper.DONORS:
            if donor in REUSED:
                continue
            guard(frozen); guard(reuse_guards)
            base = output/donor; base.mkdir(); (base/'raw').mkdir()
            bindings, filters, acquired, aligns = {}, {}, [], []
            helper.bound(output, base, donor_caps[donor])
            for row in sorted((r for r in rows if r['donor']==donor), key=lambda r:r['run_accession']):
                phase = 'acquisition'
                fastq, aligned, filtered, align, to_bam, prefilter = commands(row, base, bowtie, index)
                acquired.append(acquire_with_retries(row, fastq, output, base, donor_caps[donor]))
                run_id = row['run_accession']
                write_json(base/(run_id+'_acquisition.json'), acquired[-1])
                phase = 'alignment'
                helper.run([align, to_bam], base/(run_id+'_alignment.log'), output, base, donor_caps[donor])
                phase = 'primary_filter'
                helper.run([prefilter], base/(run_id+'_filter.log'), output, base, donor_caps[donor])
                helper.run([[str(helper.ENV/'samtools'), 'quickcheck', str(aligned), str(filtered)]],
                           base/(run_id+'_quickcheck.log'), output, base, donor_caps[donor])
                bindings[str(filtered)] = {k:row[k] for k in
                    ('run_accession', 'donor', 'GSM', 'experiment_accession', 'source_biosample')}
                filters[str(filtered)] = prefilter
                aligns.append(dict(alignment=align, to_bam=to_bam, prefilter=prefilter))
            phase = 'fixed_region_counting'
            count_path = base/'native_featurecounts.tsv'
            command = count_command(base, bindings)
            helper.run([command], base/'featurecounts.log', output, base, donor_caps[donor])
            guard(frozen); guard(reuse_guards)
            receipt = base/'summary.json'
            write_json(receipt, dict(donor=donor, h3_reads_processed=True, count_table=str(count_path),
                       count_command=command, columns=bindings, prefilter_commands=filters,
                       alignment_commands=aligns, count_sha256=sha(count_path),
                       count_summary_sha256=sha(Path(str(count_path)+'.summary')), native_read_files=acquired,
                       saved_output_bytes=native_output_bytes(base), predictions_refitted=False,
                       count_axis_assembly_or_target_evaluation_run=False))
            receipts[donor] = receipt
            completed.append(donor)
            print(json.dumps(dict(donor=donor, completed=len(completed), total=17)), flush=True)
        guard(frozen); guard(reuse_guards)
        require(tuple(completed)==helper.DONORS and set(receipts)==set(helper.DONORS), 'Require all fixed seventeen')
        require(combined_output_bytes(output)<=helper.STOP, 'Combined storage stop exceeded')
        write_json(output/'summary.json', dict(completed=helper.DONORS, all17_reads_processed=True, biological_n=17,
                   donor_receipt_paths={d:str(receipts[d]) for d in helper.DONORS}, reused_donors=list(REUSED),
                   donor_summary_sha256={d:sha(receipts[d]) for d in helper.DONORS},
                   saved_output_bytes=combined_output_bytes(output), new_saved_output_bytes=native_output_bytes(output),
                   original_saved_output_bytes=native_output_bytes(ORIGINAL),
                   protocol_sha256=sha(output/'protocol.json'), predictive_accuracy_evaluated=False,
                   predictions_refitted=False, original_processing_source_sha256=SOURCE_SHA,
                   recovery_source_sha256=args.expected_source_sha256, count_payloads_parsed_or_copied=False))
    except BaseException:
        write_json(output/'failure.json', dict(donor=donor, phase=phase, completed=completed,
                   partials_preserved=True, original_outputs_untouched=True, traceback=traceback.format_exc()))
        raise


if __name__=='__main__':
    main()
