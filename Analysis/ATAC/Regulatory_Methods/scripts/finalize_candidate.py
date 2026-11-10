#!/usr/bin/env python3
"""Assemble existing candidate evidence; never select methods or change authorities."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def table(path):
    if not path.is_file() or path.stat().st_size < 2:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def inspect(root, kind):
    if not root.exists():
        return 'not_run', 'directory absent'
    if any(root.glob('INVALIDATED*.json')):
        return 'failed', 'invalidated; excluded'
    marker = root / ('run_summary.json' if kind == 'da' else 'completion.json')
    if not marker.is_file():
        return 'partial', 'completion summary absent'
    try:
        payload = json.loads(marker.read_text())
    except (ValueError, OSError):
        return 'failed', 'unreadable completion summary'
    if kind != 'da' and payload.get('status') != 'completed':
        return 'partial', 'completion status ' + str(payload.get('status'))
    required = ['cohort_summary.tsv', 'baseline_reproduction.tsv'] if kind == 'da' else ['comparisons.tsv', 'run_settings.json']
    if any(not (root / name).is_file() for name in required):
        return 'partial', 'required evidence table absent'
    return 'completed', 'candidate evidence available'


def validate_motif_annotation(annotation, networks, motif_root):
    """Check a completed annotation against the exact six current input contexts."""
    annotation = annotation.resolve()
    if not annotation.is_dir() or any(annotation.glob('INVALIDATED*')):
        raise ValueError('Missing or invalidated precomputed motif annotation')
    completion = json.loads((annotation / 'completion.json').read_text())
    if completion.get('status') != 'completed' or completion.get('exact_ID_only') is not True:
        raise ValueError('Precomputed motif annotation is incomplete or lacks exact-ID matching')
    if completion.get('motif_root') != str(motif_root.resolve()):
        raise ValueError('Precomputed annotation uses a different frozen motif root')
    database_records = table(motif_root / 'jaspar_database.tsv')
    if len(database_records) != 1 or completion.get('jaspar_source') != database_records[0]:
        raise ValueError('Precomputed annotation uses a different frozen JASPAR source')
    record = database_records[0]
    digest = hashlib.sha256()
    with Path(record['database_path']).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    if digest.hexdigest() != record['database_sha256']:
        raise ValueError('Frozen JASPAR database SHA256 differs')
    roots = {str(path.resolve()): path for path in networks}
    sources = completion.get('sources', [])
    if len(roots) != 6 or len(sources) != 6 or {s['network_root'] for s in sources} != set(roots):
        raise ValueError('Precomputed annotation must cover the exact six completed network roots')
    expected = {}
    for source in sources:
        network = roots[source['network_root']]
        if inspect(network, 'network')[0] != 'completed' or any(network.glob('INVALIDATED*')):
            raise ValueError('Precomputed annotation network is incomplete or invalidated')
        current = json.loads((network / 'completion.json').read_text())
        scope = 'technical_pilot' if current.get('pilot') is True else 'development_candidate' if current.get('pilot') is False else 'scope_not_declared'
        if source.get('completion') != current or source.get('source_scope') != scope or source.get('source_interpretation') != current.get('interpretation', 'not_declared'):
            raise ValueError('Precomputed annotation network completion/source identity differs')
        comparisons = {(row['cohort'], row['lineage']): row for row in table(network / 'comparisons.tsv')}
        contexts = [d for d in network.iterdir() if d.is_dir() and (d / 'links.tsv').is_file()]
        if len(contexts) != 1:
            raise ValueError('Require one link context per completed full network root')
        for directory in contexts:
            if any(directory.glob('INVALIDATED*')):
                raise ValueError('Invalidated link context')
            cohort, lineage = directory.name.split('_', 1)
            with (directory / 'links.tsv').open() as handle:
                header = next(csv.reader(handle, delimiter='\t'))
            if not {'outer_fold','gene_id','peak_id','method','rank','score','score_unit','evidence_state'}.issubset(header):
                raise ValueError('Incomplete source link schema')
            meta = comparisons.get((cohort, lineage), {})
            expected[(source['network_root'], cohort, lineage)] = dict(source_scope=scope,
                source_interpretation=current.get('interpretation', 'not_declared'),
                lineage_state=meta.get('lineage_state', 'source_label'),
                source_lineage=meta.get('source_lineage', lineage), endpoint_recorded=str('endpoint' in header))
    coverage = table(annotation / 'coverage.tsv')
    if len(coverage) != 6 or len(expected) != 6:
        raise ValueError('Precomputed motif coverage is missing or repeated')
    if {(key[1], key[2]) for key in expected} != {(c, l) for c in ('GSE244832', 'GSE281367', 'GSE296875') for l in ('hepatocyte', 'stellate')}:
        raise ValueError('Precomputed motif annotation does not cover all six planned contexts')
    seen = set()
    for row in coverage:
        key = (row['network_root'], row['cohort'], row['lineage'])
        if key in seen or key not in expected or any(row.get(k) != v for k, v in expected[key].items()):
            raise ValueError('Precomputed motif coverage source/context differs')
        seen.add(key)
        directory = motif_root / 'motif_matches' / row['cohort'] / row['lineage']
        available = row['cohort'] != 'GSE296875' and all((directory / name).is_file() for name in ('matches.mtx','peaks.tsv','motifs.tsv'))
        reason = '' if available else 'whole_native_peak_universe_motifscan_not_run_125_peak_pilot_not_global' if row['cohort'] == 'GSE296875' else 'native_peak_motifscan_not_run'
        if row['status'] != ('annotated' if available else 'untestable') or row['reason'] != reason or row['peak_matching'] != 'exact_ID_only_no_overlap_borrowing':
            raise ValueError('Incorrect precomputed motif coverage evidence state')
        if int(row['n_link_records']) < 1 or int(row['n_family_gene_rows']) < 0 or (not available and int(row['n_family_gene_rows']) != 0):
            raise ValueError('Incomplete precomputed motif coverage counts')
    with gzip.open(annotation / 'motif_family_gene.tsv.gz', 'rt') as handle:
        header = next(csv.reader(handle, delimiter='\t'))
    if not {'network_root','cohort','lineage','gene_id','motif_family','score_unit','evidence_state','tf_activity_inferred','tf_binding_inferred','causal_regulation_inferred','jaspar_database_sha256'}.issubset(header):
        raise ValueError('Missing precomputed motif-family table/schema')
    return dict(annotation=str(annotation), completion=completion,
        completion_sha256=hashlib.sha256((annotation / 'completion.json').read_bytes()).hexdigest(),
        coverage_sha256=hashlib.sha256((annotation / 'coverage.tsv').read_bytes()).hexdigest())


def assemble(execution, output, source, rscript, render=True, motifs=True, motif_annotation=None):
    output.mkdir(parents=True, exist_ok=False)
    motif_root = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/chromvar')
    snapshot = output / 'source_snapshot'
    snapshot.mkdir()
    hashes = {}
    for name in ('finalize_candidate.py', 'summarize_candidate.py', 'render_results.R',
                 'annotate_motif_families.py', 'annotate_motif_context.py'):
        path = source / name
        if path.is_file():
            content = path.read_bytes()
            (snapshot / name).write_bytes(content)
            hashes[str(path.resolve())] = hashlib.sha256(content).hexdigest()
    (output / 'source_sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')
    source = snapshot
    records, networks, das = [], [], {}
    for name in ('full', 'cell_downsample'):
        path = execution / 'da_full' / name
        status, reason = inspect(path, 'da')
        records.append(dict(component='da_' + name, path=str(path.resolve()), status=status, reason=reason))
        if status == 'completed':
            das[name] = path
    for cohort in ('GSE244832', 'GSE281367', 'GSE296875'):
        for lineage in ('hepatocyte', 'stellate'):
            path = execution / 'networks_full' / f'{cohort}_{lineage}' / 'analysis'
            status, reason = inspect(path, 'network')
            records.append(dict(component=f'{cohort}_{lineage}', path=str(path.resolve()), status=status, reason=reason))
            if status == 'completed':
                networks.append(path)
    motif_receipt = None
    if motif_annotation is not None:
        if not motifs:
            raise ValueError('--motif-annotation cannot be combined with --skip-motifs')
        motif_receipt = validate_motif_annotation(motif_annotation, networks, motif_root)
    def run(label, command):
        with (output / (label + '.log')).open('w') as handle:
            result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        records.append(dict(component=label, path=str(output.resolve()), status='completed' if result.returncode == 0 else 'failed', reason=f'exit {result.returncode}; see {label}.log'))
        return result.returncode == 0
    # No substitute DA is fabricated when the full DA run is unfinished.
    if 'full' in das:
        command = [sys.executable, str(source / 'summarize_candidate.py'), '--da', str(das['full']), '--output', str(output / 'pooled')]
        for path in networks:
            command.extend(['--network', str(path)])
        run('fixed_family_summary', command)
    else:
        records.append(dict(component='fixed_family_summary', path='', status='not_run', reason='full DA incomplete; no substitute DA'))
    if render:
        for name, path in das.items():
            run('figures_da_' + name, [str(rscript), str(source / 'render_results.R'), '--da', str(path), '--output', str(output / ('figures_da_' + name))])
        for path in networks:
            name = path.parent.name
            run('figures_' + name, [str(rscript), str(source / 'render_results.R'), '--network', str(path), '--output', str(output / ('figures_' + name))])
    if motif_annotation is not None:
        receipt = motif_receipt
        (output / 'motif_family_compatibility').symlink_to(Path(receipt['annotation']), target_is_directory=True)
        (output / 'motif_annotation_reference.json').write_text(json.dumps(receipt, indent=2) + '\n')
        records.append(dict(component='motif_family_compatibility', path=receipt['annotation'],
            status='completed', reason='validated read-only precomputed annotation; see motif_annotation_reference.json'))
    elif motifs and networks:
        annotation = source / 'annotate_motif_families.py'
        if annotation.exists():
            command = [sys.executable, str(annotation), '--motif-root', str(motif_root), '--output', str(output / 'motif_family_compatibility')]
            for path in networks:
                command.extend(['--network', str(path)])
            run('motif_family_compatibility', command)
        else:
            records.append(dict(component='motif_family_compatibility', path='', status='not_run', reason='annotation script unavailable'))
    lines = ['Candidate snATAC comparison. No scientific adoption or automatic method selection.',
             'Completed means available analytical output, not biological validation.', '']
    for row in records:
        lines.append(f"{row['component']}: {row['status']}; {row['reason']}")
    for name, path in das.items():
        lines.append(f"DA {name}: {len(table(path/'cohort_summary.tsv'))} cohort/lineage/representation/filter/method rows; "
                     f"{len(table(path/'null_diagnostics.tsv'))} null diagnostic rows; {len(table(path/'lodo_summary.tsv'))} donor-omission rows.")
        for label, filename in (('Null diagnostics', 'null_diagnostics.tsv'), ('Donor omission', 'lodo_summary.tsv')):
            diagnostics = table(path / filename)
            statuses = {}
            for row in diagnostics:
                status = row.get('status') or 'status_unrecorded'
                statuses[status] = statuses.get(status, 0) + 1
            lines.append(f"DA {name} {label}: " + (json.dumps(statuses, sort_keys=True) if diagnostics else 'unavailable; no diagnostic rows'))
    for path in networks:
        for directory in path.iterdir():
            if not directory.is_dir():
                continue
            coverage = table(directory / 'coverage.tsv')
            if coverage:
                settings = json.loads((path / 'run_settings.json').read_text())
                counts = {}
                for row in coverage:
                    counts[row['status']] = counts.get(row['status'], 0) + 1
                lines.append(f"Context {path.parent.name}/{directory.name}: scope={settings.get('analysis_scope','unrecorded')}; "
                             f"max_genes={settings.get('max_genes','unrecorded')}; gene-fold coverage={counts}.")
    pooled = table(output / 'pooled' / 'network_comparisons.tsv')
    if pooled:
        evaluated = [r for r in pooled if r['status'] == 'evaluated']
        supported = [r for r in evaluated if r.get('predictive_improvement','').lower() == 'true']
        lines.append(f"Fixed 48-comparison Holm family: {len(evaluated)} evaluated; {len(supported)} satisfy positive donor-bootstrap interval and Holm support in nontechnical analyses.")
        for row in supported:
            lines.append('Conditional support: ' + ' | '.join(f'{key}={row.get(key, "")}' for key in ('cohort','lineage','endpoint','method','baseline','improvement','ci_low','ci_high','holm_pvalue','n_donors','n_common_genes')))
    lines.extend(['', 'Prediction uncertainty is conditional on fixed donor cross-validation predictions; donors are the inference unit.',
                  'DA nulls are fitted-intercept negative-binomial diagnostics on a selected subset of at most 2,000 fixed peaks; they do not establish biological truth or full-universe BH calibration.',
                  'Binary accessible-nucleus nulls also assume negative-binomial counts; they do not establish calibration under the true bounded binary-aggregation distribution.',
                  'More significant peaks do not establish superiority. Current voom remains the manuscript reference; assess effects, uncertainty, cohort agreement, donor influence and calibration jointly.',
                  'Method scores, accessibility links and motif-family compatibility do not establish TF activity or causal regulation.',
                  'FNIH public donor matrices and independence remain unqualified here; no external MASH replication claim is added.',
                  'GSE244832 RNA identity/repeated-run resolution remains required; its separate libraries are not same-nucleus validation.',
                  'Native TRIPOD/LINGER work is separate bounded development, not a completed full multimodal comparison.',
                  'GSE296875 mesenchymal source labels do not establish stellate-specific validation.',
                  'Current manuscript figures, frozen programs and Gene Catalog evidence are unchanged.'])
    (output / 'REPORT.txt').write_text('\n'.join(lines) + '\n')
    with (output / 'components.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=['component','path','status','reason'], delimiter='\t'); writer.writeheader(); writer.writerows(records)
    (output / 'assembly.json').write_text(json.dumps({'candidate_only':True,'scientific_adoption':False,'execution_root':str(execution.resolve()),'components':records},indent=2)+'\n')
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execution-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source-root',type=Path,default=Path(os.environ.get('REGULATORY_SOURCE_ROOT',Path(__file__).parent)))
    parser.add_argument('--rscript',type=Path,default=Path('/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript'))
    parser.add_argument('--skip-render',action='store_true'); parser.add_argument('--skip-motifs',action='store_true')
    parser.add_argument('--motif-annotation', type=Path, help='Reuse an exact-source completed motif annotation read-only')
    args=parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Execute on an allocated compute node')
    records = assemble(args.execution_root,args.output,args.source_root,args.rscript,not args.skip_render,not args.skip_motifs,args.motif_annotation)
    if any(row['status'] == 'failed' for row in records):
        raise RuntimeError('Candidate assembly contains failed components; inspect components.tsv')


if __name__ == '__main__':
    main()
