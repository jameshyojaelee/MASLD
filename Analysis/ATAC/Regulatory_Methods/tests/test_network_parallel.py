"""Bounded fork execution preserves annotation order and reports worker failures."""
import importlib.util
import os
from pathlib import Path
import time

import pytest

SPEC = importlib.util.spec_from_file_location('unregistered_parallel_core',
    Path(__file__).resolve().parents[1]/'scripts/network_compare.py')
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_worker_counts_cannot_exceed_allocation(monkeypatch):
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '2')
    assert module.worker_count(2, 'test') == 2
    for value in [0, -1, 3]:
        with pytest.raises(ValueError):
            module.worker_count(value, 'test')


def test_importlib_worker_order_and_readonly_context(monkeypatch):
    assert os.environ.get('SLURM_JOB_ID')
    def calculate(job, context):
        if job[0] == 1:
            time.sleep(.02)
        with pytest.raises(TypeError):
            context['value'] = -1
        return job[0], context['value']
    monkeypatch.setattr(module, '_fit_gene_readout', calculate)
    jobs = [(i, {'gene_id': str(i)}) for i in range(1, 15)]
    original = {'value': 42}
    assert list(module._ordered_gene_readouts(jobs, original, 2)) == [(i, 42) for i in range(1, 15)]
    assert original == {'value': 42}


def test_worker_failure_identifies_gene_and_reaches_parent(monkeypatch):
    assert os.environ.get('SLURM_JOB_ID')
    def calculate(job, context):
        raise ValueError('fixture failure')
    monkeypatch.setattr(module, '_fit_gene_readout', calculate)
    with pytest.raises(RuntimeError, match='gene bad_gene'):
        list(module._ordered_gene_readouts([(1, {'gene_id': 'bad_gene'})], {}, 2))
