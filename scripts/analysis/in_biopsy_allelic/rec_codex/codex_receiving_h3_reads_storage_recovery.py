"""Fixed H3 processing with an explicit B36 operational retention amendment.

Both failed jobs stay in place. Count bytes are hashed only, never parsed.
No target/prediction values, model fit, calibration, selection or accuracy.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ORIGINAL = REC/'codex_receiving_h3_reads_21998882'
FAILED_RECOVERY = REC/'codex_receiving_h3_reads_recovery_22023804'
PRIOR_SOURCE = HERE/'codex_receiving_h3_reads_recovery.py'
PRIOR_SOURCE_SHA = '733f2c793206fcb122250d3f39aa2c19a9705d926858ed32bd3aa5f36180f01d'
ORIGINAL_SOURCE_SHA = 'f641ea6485515b7f986f3ce279f449c5e6b85e99f676c3a5ac959fc93bf8b509'
SOURCE_SHA = ORIGINAL_SOURCE_SHA
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
FAILED_GUARDS = {
    FAILED_RECOVERY/'protocol.json':'34e6c6cd524986d9965f1beb870b0fd8a54404837329dfa988f5fabafdd95e56',
    FAILED_RECOVERY/'failure.json':'8773a9f5a1e41d7708883fc72e8d302e716e108d65e37eeb365fd79dee5e1ab4',
    FAILED_RECOVERY/'codex_receiving_h3_reads_recovery.py':PRIOR_SOURCE_SHA,
}
B36_CAP = 40_000_000_000
OLD_B36_CAP = 22_047_488_435
GLOBAL_STOP = 230_000_000_000
GLOBAL_ACQUISITION_CAP = 250_000_000_000


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source-sha256',required=True)
    args = parser.parse_args()
    require(re.fullmatch('[0-9a-f]{64}',args.expected_source_sha256) is not None
            and sha(Path(__file__))==args.expected_source_sha256,'Submitted source differs')
    require(os.environ.get('SLURM_JOB_ID','').isdigit()
            and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=8,'Require eight allocated CPUs')
    require(sha(PRIOR_SOURCE)==PRIOR_SOURCE_SHA
            and sha(HERE/'codex_receiving_h3_reads.py')==ORIGINAL_SOURCE_SHA,
            'Original/prior processing helpers changed')
    # Import only the prospectively hash-verified fixed processing definitions.
    import codex_receiving_h3_reads_recovery as prior
    helper = prior.helper
    require(prior.ORIGINAL==ORIGINAL and prior.SOURCE_SHA==ORIGINAL_SOURCE_SHA
            and prior.ORIGINAL_GUARDS==ORIGINAL_GUARDS and prior.REUSED==REUSED
            and helper.STOP==GLOBAL_STOP and helper.CAP==GLOBAL_ACQUISITION_CAP,
            'Prior roots, source or unchanged global allowances differ')
    helper.guard(FAILED_GUARDS)
    failed_protocol = helper.read_json(FAILED_RECOVERY/'protocol.json')
    failure = helper.read_json(FAILED_RECOVERY/'failure.json')
    require(not (FAILED_RECOVERY/'summary.json').exists() and failure['donor']=='B36'
            and failure['phase']=='alignment' and tuple(failure['completed'])==tuple(prior.REUSED)
            and failure['partials_preserved'] is True and failure['original_outputs_untouched'] is True,
            'Prior failed storage boundary changed')
    require(failed_protocol['recovery_source_sha256']==PRIOR_SOURCE_SHA
            and failed_protocol['original_processing_directory']==str(ORIGINAL)
            and failed_protocol['donor_caps']['B36']==OLD_B36_CAP
            and failed_protocol['sampled_stop_bytes']==GLOBAL_STOP
            and failed_protocol['hard_cap_bytes']==GLOBAL_ACQUISITION_CAP
            and failed_protocol['count_payloads_parsed_or_copied'] is False
            and failed_protocol['storage_scope']=='all original retained files plus all new files; per-donor also combines old/new',
            'Prior operational storage protocol differs')
    frozen, rows, index = helper.prerequisites()
    bowtie = shutil.which('bowtie2')
    require(bowtie,'Load Bowtie2/2.5.4-linux-x86_64')
    versions = {}
    for name, cmd in [('bowtie2',[bowtie,'--version']),
                      ('samtools',[str(helper.ENV/'samtools'),'--version']),
                      ('featureCounts',[str(helper.ENV/'featureCounts'),'-v'])]:
        p = subprocess.run(cmd,capture_output=True,text=True,check=True)
        versions[name] = p.stdout+p.stderr
    require(versions==failed_protocol['versions'],'Tool versions changed from both prior runs')
    frozen.update({Path(bowtie):sha(Path(bowtie)),helper.ENV/'samtools':sha(helper.ENV/'samtools'),
                   helper.ENV/'featureCounts':sha(helper.ENV/'featureCounts')})
    reuse_guards, receipts = prior.admit_reused(frozen,rows,versions,bowtie,index)
    frozen.update(FAILED_GUARDS)
    frozen.update({PRIOR_SOURCE:PRIOR_SOURCE_SHA,HERE/'codex_receiving_h3_reads.py':ORIGINAL_SOURCE_SHA,
                   HERE/'run_codex_receiving_h3_reads.sbatch':sha(HERE/'run_codex_receiving_h3_reads.sbatch'),
                   Path(__file__):args.expected_source_sha256,
                   HERE/'run_codex_receiving_h3_reads_storage_recovery.sbatch':
                       sha(HERE/'run_codex_receiving_h3_reads_storage_recovery.sbatch')})
    output = REC/('codex_receiving_h3_reads_storage_recovery_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    native_output_bytes = helper.output_bytes
    retained_roots = (ORIGINAL,FAILED_RECOVERY)

    def combined_output_bytes(path):
        size = native_output_bytes(path)
        if path==output:
            size += sum(native_output_bytes(root) for root in retained_roots)
        elif path.parent==output and path.name in helper.DONORS:
            size += sum(native_output_bytes(root/path.name) for root in retained_roots
                        if (root/path.name).exists())
        return size

    helper.output_bytes = combined_output_bytes
    donor_caps = {d:5*sum(int(r['fastq_bytes']) for r in rows if r['donor']==d)+2_000_000_000
                  for d in helper.DONORS}
    require(donor_caps==failed_protocol['donor_caps'] and donor_caps['B36']==OLD_B36_CAP,
            'Prior donor storage allowance formula changed')
    donor_caps['B36'] = B36_CAP
    amendment = dict(donor='B36',previous_cap_bytes=OLD_B36_CAP,amended_cap_bytes=B36_CAP,
                     all_other_donor_caps_unchanged=True,global_sampled_stop_bytes=GLOBAL_STOP,
                     global_acquisition_cap_bytes=GLOBAL_ACQUISITION_CAP,campaign_allowance_bytes=1_500_000_000_000,
                     reason='Preserve both failed jobs while repeating fixed B36 processing; operational retention only',
                     estimated_extra_B36_bytes=13_000_000_000,estimated_combined_B36_bytes=35_000_000_000,
                     estimate_label='Estimate, not measured required storage; no automatic cap increase')
    protocol = dict(helper.read_json(ORIGINAL/'protocol.json'))
    protocol.update(input_sha256={str(p):h for p,h in frozen.items()},
                    reused_input_sha256={str(p):h for p,h in reuse_guards.items()},
                    donor_caps=donor_caps,operational_cap_amendment=amendment,
                    original_processing_directory=str(ORIGINAL),prior_failed_recovery_directory=str(FAILED_RECOVERY),
                    retained_processing_directories=[str(root) for root in retained_roots],
                    reused_donors=list(prior.REUSED),reacquired_donors=[d for d in helper.DONORS if d not in prior.REUSED],
                    recovery_source_sha256=args.expected_source_sha256,prior_recovery_source_sha256=PRIOR_SOURCE_SHA,
                    download_max_attempts=3,retry_delay_seconds=[10,20],
                    retry_partial_policy='Preserve each failed transport attempt; restart same full URL without Range',
                    count_payloads_reused_in_place=True,count_payloads_parsed_or_copied=False,
                    storage_scope='All files in both retained prior directories plus new output, globally and per donor',
                    storage_failure_policy='Stop and preserve all partials; no automatic storage retry or cap increase')
    helper.write_json(output/'protocol.json',protocol)
    (output/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    helper.write_json(output/'reused_receipts.json',dict(
        donor_receipt_paths={d:str(p) for d,p in receipts.items()},donor_summary_sha256=prior.REUSED,
        original_source_sha256=ORIGINAL_SOURCE_SHA))
    completed = list(prior.REUSED)
    donor, phase = None, 'initialization'
    try:
        for donor in helper.DONORS:
            if donor in prior.REUSED:
                continue
            helper.guard(frozen);helper.guard(reuse_guards)
            base = output/donor;base.mkdir();(base/'raw').mkdir()
            bindings, filters, acquired, aligns = {}, {}, [], []
            helper.bound(output,base,donor_caps[donor])
            for row in sorted((r for r in rows if r['donor']==donor),key=lambda r:r['run_accession']):
                phase = 'acquisition'
                fastq, aligned, filtered, align, to_bam, prefilter = prior.commands(row,base,bowtie,index)
                acquired.append(prior.acquire_with_retries(row,fastq,output,base,donor_caps[donor]))
                run_id = row['run_accession']
                helper.write_json(base/(run_id+'_acquisition.json'),acquired[-1])
                phase = 'alignment'
                helper.run([align,to_bam],base/(run_id+'_alignment.log'),output,base,donor_caps[donor])
                phase = 'primary_filter'
                helper.run([prefilter],base/(run_id+'_filter.log'),output,base,donor_caps[donor])
                helper.run([[str(helper.ENV/'samtools'),'quickcheck',str(aligned),str(filtered)]],
                           base/(run_id+'_quickcheck.log'),output,base,donor_caps[donor])
                bindings[str(filtered)] = {k:row[k] for k in
                    ('run_accession','donor','GSM','experiment_accession','source_biosample')}
                filters[str(filtered)] = prefilter
                aligns.append(dict(alignment=align,to_bam=to_bam,prefilter=prefilter))
            phase = 'fixed_region_counting'
            count_path = base/'native_featurecounts.tsv'
            command = prior.count_command(base,bindings)
            helper.run([command],base/'featurecounts.log',output,base,donor_caps[donor])
            helper.guard(frozen);helper.guard(reuse_guards)
            receipt = base/'summary.json'
            helper.write_json(receipt,dict(donor=donor,h3_reads_processed=True,count_table=str(count_path),
                count_command=command,columns=bindings,prefilter_commands=filters,alignment_commands=aligns,
                count_sha256=sha(count_path),count_summary_sha256=sha(Path(str(count_path)+'.summary')),
                native_read_files=acquired,saved_output_bytes=native_output_bytes(base),predictions_refitted=False,
                count_axis_assembly_or_target_evaluation_run=False))
            helper.bound(output,base,donor_caps[donor])
            receipts[donor] = receipt;completed.append(donor)
            print(json.dumps(dict(donor=donor,completed=len(completed),total=17)),flush=True)
        helper.guard(frozen);helper.guard(reuse_guards)
        require(tuple(completed)==helper.DONORS and set(receipts)==set(helper.DONORS),'Require all fixed seventeen')
        require(combined_output_bytes(output)<=GLOBAL_STOP,'Combined sampled storage stop exceeded')
        helper.write_json(output/'summary.json',dict(completed=helper.DONORS,all17_reads_processed=True,
            biological_n=17,donor_receipt_paths={d:str(receipts[d]) for d in helper.DONORS},
            donor_summary_sha256={d:sha(receipts[d]) for d in helper.DONORS},reused_donors=list(prior.REUSED),
            saved_output_bytes=combined_output_bytes(output),new_saved_output_bytes=native_output_bytes(output),
            retained_saved_output_bytes_by_directory={str(root):native_output_bytes(root) for root in retained_roots},
            protocol_sha256=sha(output/'protocol.json'),operational_cap_amendment=amendment,
            predictive_accuracy_evaluated=False,predictions_refitted=False,count_payloads_parsed_or_copied=False,
            original_processing_source_sha256=ORIGINAL_SOURCE_SHA,prior_recovery_source_sha256=PRIOR_SOURCE_SHA,
            recovery_source_sha256=args.expected_source_sha256,all_prior_files_preserved=True))
        require(combined_output_bytes(output)<=GLOBAL_STOP,'Final combined sampled storage stop exceeded')
    except BaseException:
        helper.write_json(output/'failure.json',dict(donor=donor,phase=phase,completed=completed,
            partials_preserved=True,all_prior_outputs_untouched=True,operational_cap_amendment=amendment,
            traceback=traceback.format_exc()))
        raise


if __name__=='__main__':
    main()
