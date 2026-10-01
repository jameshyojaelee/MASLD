#!/usr/bin/env python3
"""One public ATAC run, first 100 spots: format metadata only, no alignment."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import selectors
import signal
import subprocess
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
META = REC / 'codex_gse223843_raw_read_metadata_20261001'
SRA = Path('/nfs/sw/easybuild/software/sratoolkit/3.2.1')
ARC = Path('/nfs/sw/easybuild/software/cellranger-arc-2.0.2')
EXPORT = SRA / 'bin/fastq-dump-orig.3.2.1'
VDB = SRA / 'bin/vdb-config.3.2.1'
WHITELIST = ARC / 'lib/python/atac/barcodes/737K-arc-v1.txt.gz'
LAUNCHER = HERE / 'run_codex_gse223843_read_layout_pilot.sbatch'
RUN = 'SRR23252452'
URL = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23252452/SRR23252452'
SRA_BYTES = 5191054225
SRA_MD5 = '3cdda9ebd3fcd6c156307fab5b98bbe6'
SAVED_CAP = NETWORK_CAP = 6000000000
FASTQ_CAP, META_CAP = 1024**2, 32 * 1024**2
SPOTS = 100
LENGTHS = {1: 8, 2: 50, 3: 24, 4: 49}
GUARDS = {
    META / 'identity_read_layout_file_metadata.json': 'daf440a0425678ccc597f1fe0e2e898949d6cd80266fb5e959c61ecf4fb01b34',
    META / 'ENA_ATAC_file_metadata.json': 'd4b48f024808f7ba5cefe963983abadb2c8918ad7be32e641ca88967311ca667',
    META / 'normalized_SRA_HEAD_metadata.json': '418fcffed5263f8e86e803111bbd7df0bcb9806e1ee3b5e47dfe36be2a255a31',
    REC / 'codex_gse223843_identity_metadata/identity_metadata.json': '172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018',
    ROOT / 'docs/technical/agent_exchange/2026-10-01_codex_gse223843_raw_fragment_feasibility.md': '6e45e6bc6a5af2b63354232186ac9e5810f8ea8ea46e31284a5830a0b9929bf8',
    EXPORT: '9f89c4a6af231c09c180b40de5b774e397a659dd79c6ca71f9b2cfbf9a7915ba',
    VDB: '9a1de9a9092fc6214cab3178df2b705efd32b205c29892831785420b12d8178c',
    SRA / 'bin/ncbi/default.kfg': 'cdd345114ddda0a1784db04938dbac2d66346a290ef8c9d383d58892760035c3',
    ARC / 'bin/cellranger-arc': '61e9fd4936c3c5bf24ca78d7f2d777ea6bafa2ec9241c265fe8e1801ca412f72',
    ARC / 'lib/python/atac/barcodes/__init__.py': 'c896d07f2bde284d6dd36a618636e2c83d9b131a2195af627dc9eebd8c9069ba',
    ARC / 'lib/python/atac_rna/barcodes.py': '9f4e0c90eec48ac8bbadd0bd69377c8fca9b1a17cac59b216283dd50fc82f185',
    ARC / 'mro/atac_rna/sc_atac_gex_counter_cs.mro': '283bcc35bca83cd819d48a7b35d449bbed1bec2e3043f15089f6f4c50152c0ac',
    WHITELIST: '9a4c8f29a05f59b9a4d7525a2f44d1778ccfa540a8832ad9f5f84047f374cc26',
}
FLAGS = dict(raw_GEX_read=False, matrices_read=False, clinical_values_read=False,
    genotype_values_inspected=False, allelic_analysis=False, restricted_inputs_read=False,
    genomic_alignment=False, count_production=False,
    model_fit_or_score=False, orientation_selected=False, barcode_correction=False,
    barcode_translation=False, lineage_assignment=False, person_independence_verified=False,
    model_or_evaluation_admission=False, sequence_or_quality_strings_reported=False)


class PilotFailure(Exception):
    """Only fixed, payload-free reason codes may enter reports."""


def require(condition, code):
    if not condition:
        raise PilotFailure(code)


class Pilot:
    def __init__(self, output, job):
        self.output = output
        self.logs = [REC / f'1592_raw_layout_{job}.{suffix}' for suffix in ('out', 'err')]
        self.started = time.monotonic()
        self.phase = 'initial_guards'
        self.network_bytes = 0
        self.commands = []
        self.peak_saved = self.peak_fastq = self.peak_meta = 0
        self.phase_receipts = {}

    def clock(self):
        require(time.monotonic() - self.started < 3300, 'pilot_time_cap')

    def sizes(self):
        saved = fastq = metadata = 0
        for path in self.output.rglob('*'):
            require(not path.is_symlink(), 'output_symlink_forbidden')
            if not path.is_file():
                continue
            n = path.stat().st_size
            saved += n
            if path.parent == self.output / 'fastq':
                require(path.name in {f'{RUN}_{i}.fastq' for i in LENGTHS},
                        'unexpected_export_filename_contents_suppressed')
                fastq += n
            elif path.parent != self.output / 'download':
                metadata += n
        for path in self.logs:
            require(not path.is_symlink(), 'SLURM_log_symlink_forbidden')
            if path.is_file():
                n = path.stat().st_size
                saved += n; metadata += n
        self.peak_saved = max(self.peak_saved, saved)
        self.peak_fastq = max(self.peak_fastq, fastq)
        self.peak_meta = max(self.peak_meta, metadata)
        require(saved <= SAVED_CAP, 'inclusive_saved_byte_cap')
        require(fastq <= FASTQ_CAP, 'aggregate_FASTQ_byte_cap')
        require(metadata <= META_CAP, 'extra_metadata_byte_cap')
        return dict(saved_bytes=saved, fastq_bytes=fastq, extra_metadata_bytes=metadata)

    def monitor(self):
        self.clock()
        return self.sizes()

    def write_json(self, name, value):
        data = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
        sizes = self.sizes()
        require(sizes['saved_bytes'] + len(data) <= SAVED_CAP
                and sizes['extra_metadata_bytes'] + len(data) <= META_CAP,
                'receipt_byte_cap')
        with (self.output / name).open('xb') as target:
            target.write(data)

    def hashes(self, path, md5=False):
        hashes = {'sha256': hashlib.sha256()}
        if md5:
            hashes['md5'] = hashlib.md5()
        size = 0
        with path.open('rb') as source:
            for block in iter(lambda: source.read(1024**2), b''):
                self.clock()
                size += len(block)
                for h in hashes.values():
                    h.update(block)
        return dict(bytes=size, **{key: h.hexdigest() for key, h in hashes.items()})

    def guard_sources(self, script_sha, previous=None):
        values = {}
        for path, expected in {**GUARDS, Path(__file__).resolve(): script_sha}.items():
            got = self.hashes(path)
            require(got['sha256'] == expected, 'pinned_source_hash_mismatch')
            values[str(path)] = got
        values[str(LAUNCHER)] = self.hashes(LAUNCHER)
        if previous is not None:
            require(values == previous, 'source_inventory_changed')
        return values

    def run(self, argv, env, timeout, capture=False):
        """Drain/hash/discard diagnostics; never save raw export stdout/stderr."""
        self.monitor()
        stdout, stderr = hashlib.sha256(), hashlib.sha256()
        streams = {0: dict(bytes=0, hash=stdout), 1: dict(bytes=0, hash=stderr)}
        captured = bytearray()
        receipt = dict(argv=[str(x) for x in argv], stdout_contents_saved=False,
                       stderr_contents_saved=False, sequence_strings_saved=False)
        self.commands.append(receipt)

        def limits():
            # Export emits exactly four named files; each has a hard file-size cap.
            resource.setrlimit(resource.RLIMIT_FSIZE, (FASTQ_CAP // 4, FASTQ_CAP // 4))
            resource.setrlimit(resource.RLIMIT_CORE, (0, 0))

        child = subprocess.Popen([str(x) for x in argv], cwd=self.output,
            env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True, preexec_fn=limits)
        began = time.monotonic()
        sel = selectors.DefaultSelector()
        sel.register(child.stdout, selectors.EVENT_READ, 0)
        sel.register(child.stderr, selectors.EVENT_READ, 1)
        try:
            while sel.get_map() or child.poll() is None:
                self.monitor()
                require(time.monotonic() - began < timeout, 'tool_time_cap')
                for key, _ in sel.select(0.05):
                    block = os.read(key.fileobj.fileno(), 65536)
                    if not block:
                        sel.unregister(key.fileobj)
                        continue
                    stream = streams[key.data]
                    stream['bytes'] += len(block); stream['hash'].update(block)
                    require(sum(x['bytes'] for x in streams.values()) <= META_CAP,
                            'discarded_tool_diagnostic_byte_cap')
                    if capture and key.data == 0:
                        require(len(captured) + len(block) <= 128 * 1024,
                                'metadata_tool_capture_cap')
                        captured.extend(block)
            receipt['returncode'] = child.wait()
            require(receipt['returncode'] == 0, 'tool_failed_diagnostics_suppressed')
        finally:
            if child.poll() is None:
                try:
                    os.killpg(child.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    child.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(child.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    child.wait(timeout=2)
            receipt['returncode'] = child.returncode
            receipt['elapsed_seconds'] = time.monotonic() - began
            for name, data in zip(('stdout', 'stderr'), streams.values()):
                receipt[name + '_bytes'] = data['bytes']
                receipt[name + '_sha256'] = data['hash'].hexdigest()
            sel.close(); child.stdout.close(); child.stderr.close()
        self.monitor()
        return bytes(captured)

    def download(self):
        class NoRedirect(urllib.request.HTTPRedirectHandler):
            def redirect_request(self, request, fp, code, msg, headers, newurl):
                raise PilotFailure('HTTP_redirect_forbidden')

        # No default proxy/fallback handler; one allowlisted HTTPS URL only.
        opener = urllib.request.build_opener(urllib.request.ProxyHandler({}), NoRedirect())
        request = urllib.request.Request(URL, headers={'Accept-Encoding': 'identity'})
        partial = self.output / 'download' / f'{RUN}.partial'
        receipt = dict(url=URL, requested_method='GET', redirects_allowed=False,
                       fallback_allowed=False, expected_bytes=SRA_BYTES,
                       expected_md5=SRA_MD5, response_body_bytes=0)
        self.phase_receipts['download'] = receipt
        began = time.monotonic()
        sha, md5 = hashlib.sha256(), hashlib.md5()
        with opener.open(request, timeout=30) as response:
            require(response.status == 200 and response.geturl() == URL, 'HTTPS_identity_status')
            length = response.headers.get('Content-Length')
            require(length is not None and length.isdigit() and int(length) == SRA_BYTES,
                    'HTTPS_content_length')
            require(response.headers.get('Content-Encoding', 'identity') == 'identity',
                    'HTTP_content_encoding')
            receipt.update(status=200, content_length=int(length))
            with partial.open('xb') as target:
                while True:
                    self.monitor()
                    require(time.monotonic() - began < 2700, 'download_time_cap')
                    remaining = SRA_BYTES - self.network_bytes
                    block = response.read(min(1024**2, remaining + 1))
                    if not block:
                        break
                    self.network_bytes += len(block)
                    receipt['response_body_bytes'] = self.network_bytes
                    require(self.network_bytes <= SRA_BYTES and self.network_bytes <= NETWORK_CAP,
                            'network_body_byte_cap_or_source_size')
                    require(self.sizes()['saved_bytes'] + len(block) <= SAVED_CAP,
                            'saved_cap_before_download_write')
                    target.write(block); sha.update(block); md5.update(block)
                target.flush(); os.fsync(target.fileno())
        require(self.network_bytes == SRA_BYTES and md5.hexdigest() == SRA_MD5,
                'download_size_or_MD5_mismatch')
        receipt.update(bytes=self.network_bytes, md5=md5.hexdigest(), sha256=sha.hexdigest(),
                       checksum_passed_before_export=True)
        local = self.output / 'download' / RUN
        require(not local.exists(), 'local_input_destination_exists')
        partial.rename(local)
        self.monitor()
        return local, dict(bytes=self.network_bytes, md5=md5.hexdigest(), sha256=sha.hexdigest())

    def isolated_config(self, env):
        config = self.output / 'isolated.mkfg'
        data = ('/repository/remote/disabled = "true"\n'
                '/repository/user/disabled = "true"\n'
                '/repository/site/disabled = "true"\n').encode()
        with config.open('xb') as target:
            target.write(data)
        env = dict(env, NCBI_SETTINGS=str(config))
        # Do not inherit alternative caller configuration search paths.
        for key in ('NCBI_HOME', 'VDB_CONFIG', 'KLIB_CONFIG', 'NCBI_VDB_CONFIG'):
            env.pop(key, None)
        guards = {}
        for repo in ('remote', 'user', 'site'):
            raw = self.run([VDB, '--output', 'x', f'/repository/{repo}/disabled'], env, 30, True)
            try:
                root = ET.fromstring(raw)
            except ET.ParseError:
                raise PilotFailure('effective_config_XML_unrecognized') from None
            # vdb-config 3.2.1 prints a queried leaf using its final node name
            # (source vdb-config.c:1594-1597), not the full repository tree.
            leaves = list(root.iter('disabled'))
            require(root.tag in ('disabled', 'VdbConfig') and len(leaves) == 1
                    and all(node.tag in ('disabled', 'VdbConfig') for node in root.iter())
                    and leaves[0].text is not None and not list(leaves[0])
                    and leaves[0].text.strip().lower() == 'true', 'effective_repository_not_disabled')
            guards[repo + '_disabled'] = True
        receipt = dict(path=str(config), **self.hashes(config), effective=guards,
                       user_configuration_written=False, alternate_environment_removed=True)
        self.phase_receipts['isolated_config'] = receipt
        return env, config, receipt

    def validate_fastq(self):
        self.monitor()
        files = sorted((self.output / 'fastq').iterdir())
        require({p.name for p in files} == {f'{RUN}_{i}.fastq' for i in LENGTHS},
                'four_expected_FASTQ_classes_missing_or_extra')
        axis = list(range(1, SPOTS + 1))
        receipts, indexes = [], None
        for index in LENGTHS:
            path = self.output / 'fastq' / f'{RUN}_{index}.fastq'
            ids, actual_indices, lengths = [], set(), Counter()
            ambiguous_bases = invalid_quality = 0
            selected = []
            with path.open('rb') as source:
                while True:
                    lines = [source.readline(513) for _ in range(4)]
                    if not lines[0]:
                        require(all(not x for x in lines), 'FASTQ_incomplete_terminal_record')
                        break
                    require(all(x.endswith(b'\n') and len(x) <= 512 for x in lines),
                            'FASTQ_line_limit_or_truncated_record')
                    header, seq, plus, qual = [x.rstrip(b'\r\n') for x in lines]
                    match = re.fullmatch(rb'@SRR23252452\.([0-9]+)/([1-4])', header)
                    require(match is not None and plus == b'+', 'FASTQ_header_or_separator_format')
                    spot, actual_index = map(int, match.groups())
                    require(actual_index == index, 'actual_defline_read_index_mismatch')
                    require(re.fullmatch(rb'[ACGTN]+', seq) is not None,
                            'FASTQ_non_DNA_sequence_suppressed')
                    require(len(seq) == len(qual) == LENGTHS[index], 'FASTQ_sequence_quality_lengths')
                    require(all(33 <= value <= 126 for value in qual), 'FASTQ_non_Phred33_printable_quality')
                    ids.append(spot); actual_indices.add(actual_index); lengths[len(seq)] += 1
                    ambiguous_bases += seq.count(b'N')
                    invalid_quality += sum(not 33 <= value <= 126 for value in qual)
                    if index == 3:
                        selected.append(seq)  # Only 100 x 24 bases held; never serialized.
                    require(len(ids) <= SPOTS, 'FASTQ_spot_count_cap')
                    self.monitor()
            require(ids == axis, 'FASTQ_actual_first100_spot_axis_not_complete')
            if indexes is None:
                indexes = ids
            require(ids == indexes, 'FASTQ_spot_order_not_synchronized')
            receipts.append(dict(file=path.name, **self.hashes(path), records=len(ids),
                observed_defline_read_indices=sorted(actual_indices),
                SRA_stored_index_zero_based=index - 1, lengths=dict(lengths),
                sequence_quality_lengths_equal=True, Phred33_printable_quality_records=len(ids),
                non_printable_quality_bytes=invalid_quality, ambiguous_base_bytes=ambiguous_bases))
            if index == 3:
                barcode_reads = selected
        return dict(files=receipts, synchronized_spot_ids=axis, actual_read_indices_verified=True,
                    source_stored_indices_to_defline_indices_offset=1,
                    expected_read_roles={1: 'I1_i7_sample_index', 2: 'R1_genomic_mate1',
                                         3: 'R2_i5_barcode_plus_spacer', 4: 'R3_genomic_mate2'},
                    roles_source_declared_not_alignment_verified=True), barcode_reads

    def validity(self, barcode_reads):
        allowed = set()
        decoded = lines = comments = 0
        with gzip.open(WHITELIST, 'rb') as source:
            while True:
                raw = source.readline(257)
                if not raw:
                    break  # Complete iteration also checks gzip CRC.
                decoded += len(raw); lines += 1
                require(decoded <= 32 * 1024**2 and len(raw) <= 256 and lines <= 1000000,
                        'canonical_whitelist_decode_cap')
                if b'#' in raw:  # Same explicit comment rule as pinned ATAC loader.
                    comments += 1
                    continue
                barcode = raw.strip()
                require(re.fullmatch(rb'[ACGT]{16}', barcode) is not None,
                        'canonical_ATAC_whitelist_schema')
                require(barcode not in allowed, 'canonical_ATAC_whitelist_duplicate')
                allowed.add(barcode)
                self.clock()
        require(allowed, 'canonical_ATAC_whitelist_empty')
        complement = bytes.maketrans(b'ACGTN', b'TGCAN')
        alternatives = {key: 0 for key in ('first16_forward', 'first16_reverse_complement',
                                         'last16_forward', 'last16_reverse_complement')}
        for seq in barcode_reads:
            hypotheses = dict(first16_forward=seq[:16], last16_forward=seq[-16:],
                first16_reverse_complement=seq[:16].translate(complement)[::-1],
                last16_reverse_complement=seq[-16:].translate(complement)[::-1])
            for key, value in hypotheses.items():
                alternatives[key] += value in allowed
        return dict(whitelist_path=str(WHITELIST), **self.hashes(WHITELIST),
            allowed_unique_entries=len(allowed), decoded_bytes=decoded, comment_lines=comments,
            gzip_CRC_passed=True, identity='ARC2.0.2_ATAC_737K-arc-v1',
            hypotheses={key: dict(valid_exact=count, spots=SPOTS, fraction=count / SPOTS)
                        for key, count in alternatives.items()}, selected_hypothesis=None,
            mismatch_correction_performed=False, translation_performed=False,
            first100_not_global_orientation_or_abundance_estimate=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    args = parser.parse_args()
    require(re.fullmatch(r'[0-9a-f]{64}', args.expected_script_sha256) is not None,
            'script_hash_argument')
    job = os.environ.get('SLURM_JOB_ID', '')
    require(re.fullmatch(r'[0-9]+', job) is not None, 'compute_allocation_required')
    output = REC / f'codex_gse223843_read_layout_pilot_{job}'
    own_output = False
    pilot = Pilot(output, job)
    try:
        output.mkdir(exist_ok=False); own_output = True
        (output / 'download').mkdir(); (output / 'fastq').mkdir()
        before = pilot.guard_sources(args.expected_script_sha256)
        with (output / 'executed_source.py').open('xb') as target:
            target.write(Path(__file__).read_bytes())
        metadata = json.loads((META / 'identity_read_layout_file_metadata.json').read_text())
        rows = [r for r in metadata['runs'] if r['run_attributes']['accession'] == RUN]
        require(len(rows) == 1, 'pinned_exact_run_missing')
        row = rows[0]
        require(row['person'] == 'Fontan1' and row['GEO'] == 'GSM6997758'
                and row['experiment'] == 'SRX19197939' and row['BioSample'] == 'SAMN32937884',
                'exact_source_person_experiment_mismatch')
        original = [f['attributes'] for f in row['file_metadata']
                    if f['attributes'].get('semantic_name') == 'SRA Normalized']
        require(len(original) == 1 and original[0]['url'] == URL
                and int(original[0]['size']) == SRA_BYTES and original[0]['md5'] == SRA_MD5,
                'normalized_file_metadata_mismatch')
        require([(int(r['index']), int(r['average']), r['stdev'])
                 for r in row['stored_read_layout']] == [(i - 1, n, '0') for i, n in LENGTHS.items()],
                'pinned_stored_read_layout_mismatch')
        env = dict(os.environ, PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1')
        pilot.phase = 'isolated_offline_configuration'
        env, config, config_before = pilot.isolated_config(env)
        versions = {}
        for name, path, expected in [('fastq_dump', EXPORT, b'3.2.1'),
                                     ('vdb_config', VDB, b'3.2.1'),
                                     ('cellranger_arc', ARC / 'bin/cellranger-arc', b'2.0.2')]:
            flag = '--version'
            raw = pilot.run([path, flag], env, 30, True)
            require(expected in raw, 'installed_tool_version_mismatch')
            versions[name] = expected.decode('ascii')
        help_raw = pilot.run([EXPORT, '--help'], env, 30, True)
        help_flags = ['--split-files', '--minSpotId', '--maxSpotId', '--defline-seq', '--offset',
                      '--defline-qual', '--skip-technical', '--disable-multithreading']
        require(all(flag.encode() in help_raw for flag in help_flags), 'installed_export_help_mismatch')
        pilot.write_json('pre_download_guards.json', dict(source_hashes=before,
            identity=dict(run=RUN, person='Fontan1', GEO='GSM6997758',
                          experiment='SRX19197939', BioSample='SAMN32937884'),
            versions=versions, verified_help_flags=help_flags,
            isolated_config=config_before, guards=FLAGS))
        pilot.phase = 'single_exact_normalized_SRA_download'
        pilot.write_json('download_started.json', dict(phase=pilot.phase, run=RUN,
            URL=URL, expected_bytes=SRA_BYTES, expected_MD5=SRA_MD5,
            source_payload_is_raw_ATAC=True, inspection_scope='first100_spot_format_only'))
        local, downloaded = pilot.download()
        pilot.phase = 'pre_export_SRA_checksum'
        require(pilot.hashes(local, True) == downloaded, 'pre_export_local_SRA_hash_mismatch')
        pilot.phase = 'first100_local_spot_export'
        argv = [EXPORT, '--split-files', '--minSpotId', '1', '--maxSpotId', str(SPOTS),
                '--disable-multithreading', '--defline-seq', '@$ac.$si/$ri',
                '--defline-qual', '+', '--offset', '33', '--ncbi_error_report', 'never',
                '--outdir', output / 'fastq', local]
        pilot.write_json('export_started.json', dict(phase=pilot.phase, argv=[str(x) for x in argv],
            source_checksum=downloaded, technical_reads_skipped=False,
            remote_user_site_repositories_disabled=True))
        pilot.run(argv, env, 300)
        pilot.phase = 'FASTQ_format_and_synchronized_identity'
        layout, barcode_reads = pilot.validate_fastq()
        pilot.phase_receipts['layout'] = layout
        pilot.phase = 'canonical_ATAC_whitelist_validity_four_unselected_hypotheses'
        validity = pilot.validity(barcode_reads)
        pilot.phase_receipts['barcode_validity'] = validity
        pilot.phase = 'final_source_SRA_and_configuration_guards'
        after = pilot.guard_sources(args.expected_script_sha256, before)
        require(pilot.hashes(local, True) == downloaded, 'post_export_local_SRA_hash_mismatch')
        require(pilot.hashes(config) == {key: config_before[key] for key in ('bytes', 'sha256')},
                'isolated_configuration_changed')
        pilot.phase = 'complete_bounded_read_layout_only'
        summary = dict(status='completed', SLURM_job_id=job, phase=pilot.phase,
            source=ROW_IDENTITY, network_body_bytes=pilot.network_bytes, source_SRA=downloaded,
            raw_ATAC_payload_acquired=True, raw_ATAC_inspection_purpose='first100_read_layout_only',
            read_layout_verified=True, Python_version=sys.version,
            phase_receipts=pilot.phase_receipts, source_hashes_before=before,
            source_hashes_after=after, commands=pilot.commands, flags=FLAGS,
            elapsed_seconds=time.monotonic() - pilot.started,
            caps=dict(network=NETWORK_CAP, saved=SAVED_CAP, FASTQ=FASTQ_CAP, extra_metadata=META_CAP),
            observed_peak_bytes=dict(saved=pilot.peak_saved, FASTQ=pilot.peak_fastq,
                                     extra_metadata=pilot.peak_meta),
            sizes_before_final_receipt=pilot.monitor())
        pilot.write_json('summary.json', summary)
        pilot.monitor()
        print(json.dumps(dict(status='completed', SLURM_job_id=job,
            metadata_receipt=str(output / 'summary.json'), sequence_contents_suppressed=True)), flush=True)
        pilot.monitor()
    except Exception as error:
        # Never stringify an arbitrary parser/tool/HTTP error or show a traceback.
        failure = dict(status='failed', phase=pilot.phase, SLURM_job_id=job,
            reason=str(error) if isinstance(error, PilotFailure) else 'unexpected_error_contents_suppressed',
            exception_type=type(error).__name__, phase_receipts=pilot.phase_receipts,
            network_body_bytes=pilot.network_bytes, commands=pilot.commands, flags=FLAGS,
            partial_files_retained=True, sequence_contents_suppressed=True,
            raw_ATAC_payload_bytes_acquired=pilot.network_bytes,
            read_layout_verified='layout' in pilot.phase_receipts)
        if own_output:
            try:
                pilot.write_json('failure.json', failure)
            except Exception:
                # Retain original failure, even if output budget prevents another receipt.
                pass
        print(json.dumps(dict(status='failed', phase=pilot.phase, SLURM_job_id=job,
                              reason=failure['reason'], sequence_contents_suppressed=True)), flush=True)
        return 1
    return 0


ROW_IDENTITY = dict(run=RUN, person='Fontan1', GEO='GSM6997758',
                    experiment='SRX19197939', BioSample='SAMN32937884')

if __name__ == '__main__':
    sys.exit(main())
