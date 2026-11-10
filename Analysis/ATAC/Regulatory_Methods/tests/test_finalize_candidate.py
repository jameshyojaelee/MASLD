#!/usr/bin/env python3
import importlib.util
import json
from pathlib import Path
import tempfile
from unittest.mock import patch
from types import SimpleNamespace
import subprocess
import sys
spec=importlib.util.spec_from_file_location('finalize',Path(__file__).resolve().parents[1]/'scripts/finalize_candidate.py')
f=importlib.util.module_from_spec(spec); spec.loader.exec_module(f)
with tempfile.TemporaryDirectory() as temp:
    root=Path(temp); execution=root/'execution'; execution.mkdir()
    bad=execution/'networks_full/GSE244832_hepatocyte/analysis'; bad.mkdir(parents=True)
    (bad/'completion.json').write_text('{"status":"completed"}')
    (bad/'INVALIDATED_test.json').write_text('{}')
    assert f.inspect(bad,'network')[0]=='failed'
    partial=execution/'da_full/full'; partial.mkdir(parents=True)
    (partial/'cohort_summary.tsv').write_text('cohort\nGSE244832\n')
    assert f.inspect(partial,'da')[0]=='partial'
    output=root/'candidate'
    f.assemble(execution,output,Path(spec.origin).parent,Path('/unused'),False,False)
    report=(output/'REPORT.txt').read_text()
    assert 'da_full: partial' in report and 'GSE244832_hepatocyte: failed' in report
    assert 'no substitute DA' in report and not (output/'pooled').exists()
    assert len(json.loads((output/'assembly.json').read_text())['components'])==9
    failed_cli_output=root/'failed_cli'
    failed_cli=subprocess.run([sys.executable,spec.origin,'--execution-root',str(execution),
        '--output',str(failed_cli_output),'--skip-render','--skip-motifs'],capture_output=True,text=True)
    assert failed_cli.returncode != 0
    assert 'failed components' in failed_cli.stderr
    assert (failed_cli_output/'components.tsv').is_file()
    try: f.assemble(execution,output,Path(spec.origin).parent,Path('/unused'),False,False)
    except FileExistsError: pass
    else: raise AssertionError('Overwrite allowed')
    (partial/'run_summary.json').write_text('{}')
    (partial/'baseline_reproduction.tsv').write_text('status\npassed\n')
    (partial/'null_diagnostics.tsv').write_text('status\ntestable\nfailed\ntestable\n')
    (partial/'lodo_summary.tsv').write_text('status\nuntestable\nnot_run\n')
    complete=root/'completed_candidate'
    f.assemble(execution,complete,Path(spec.origin).parent,Path('/unused'),False,False)
    rows=f.table(complete/'pooled/network_comparisons.tsv')
    assert len(rows)==48 and all(row['status']=='not_run' for row in rows)
    assert 'Fixed 48-comparison Holm family: 0 evaluated; 0 satisfy' in (complete/'REPORT.txt').read_text()
    report=(complete/'REPORT.txt').read_text()
    assert 'Null diagnostics: {"failed": 1, "testable": 2}' in report
    assert 'Donor omission: {"not_run": 1, "untestable": 1}' in report
    assert 'at most 2,000 fixed peaks' in report and 'true bounded binary-aggregation distribution' in report
    assert 'More significant peaks do not establish superiority' in report
    assert 'Current voom remains the manuscript reference' in report
    clean=root/'motif_execution'; network=clean/'networks_full/GSE244832_hepatocyte/analysis'
    network.mkdir(parents=True)
    (network/'completion.json').write_text('{"status":"completed"}')
    (network/'comparisons.tsv').write_text('cohort\tlineage\n')
    (network/'run_settings.json').write_text('{"analysis_scope":"genome_wide"}')
    captured=[]
    def capture(command, **kwargs):
        captured.append(command)
        return SimpleNamespace(returncode=0)
    snapshot_output=root/'motif_candidate'
    with patch.object(f.subprocess,'run',side_effect=capture):
        f.assemble(clean,snapshot_output,Path(spec.origin).parent,Path('/unused'),False,True)
    assert len(captured)==1
    command=captured[0]
    assert Path(command[1])==snapshot_output/'source_snapshot/annotate_motif_families.py'
    assert command[command.index('--motif-root')+1]=='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/chromvar'
    assert (snapshot_output/'source_snapshot/annotate_motif_context.py').is_file()
    help_run=subprocess.run([sys.executable,command[1],'--help'],capture_output=True,text=True)
    assert help_run.returncode==0 and '--motif-root' in help_run.stdout
print('PASS: partial DA, invalidation exclusion, no fabricated DA, explicit missing contexts, overwrite refusal')
print('PASS: snapshotted motif-family script imports copied helper; explicit absolute frozen motif root')
