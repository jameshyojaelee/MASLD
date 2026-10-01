"""Allocated-CPU verification of the synthetic generator; no fit or P1 test."""
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import scipy

from synthetic_p1_preascertainment import generate


def verify(sample):
    d = sample.design
    assert len(d.participants) == sample.metadata['participants']
    assert d.gene_count == 8
    assert len(sample.candidate_tags) == len(d.participants)*8*2
    assert len(sample.chosen_candidates) == len(d.participants)*8
    assert sample.generating_parameters[16] == 0.
    assert sample.metadata['tau'] == 0.
    q = sample.preselection_stage_probabilities
    assert q.shape == (len(d.participants), 3)
    assert np.isfinite(q).all() and (q > 0).all()
    np.testing.assert_allclose(q.sum(axis=1), 1., atol=1e-14, rtol=0)
    # Nonconstant conditional stage probabilities establish the designed
    # correlation without demanding significance from one realized sample.
    assert np.ptp(q[:, 0]) > .1
    for c in set(d.cohort):
        np.testing.assert_allclose(d.stage_centered[d.cohort == c].mean(), 0., atol=1e-14)
        np.testing.assert_allclose(d.z[d.cohort == c].mean(axis=0), 0., atol=1e-14)
    by_key = {}
    for row in sample.candidate_tags:
        assert 0 <= row.a <= row.n
        np.testing.assert_allclose(row.prior, [.42, .49, .09], atol=1e-14, rtol=0)
        by_key.setdefault((row.individual, row.gene), []).append(row)
    for row in sample.chosen_candidates:
        tags = by_key[row.individual, row.gene]
        winner = tags[int(np.argmax([r.n for r in tags]))]
        assert row == winner  # includes lower-position choice when depths tie
    selected = {(row.individual, row.gene):row for row in d.rows}
    assert len(selected) == len(d.rows)
    for row in sample.chosen_candidates:
        assert ((row.individual, row.gene) in selected) == (row.a in row.allowed)
    for row in d.rows:
        assert row.a in row.allowed
    assert sample.metadata['zero_candidate_participants'] >= 6
    assert sample.metadata['zero_selected_participants'] >= sample.metadata['zero_candidate_participants']
    assert sum(sample.metadata['candidate_genotype_class_counts_H_R_A']) == len(d.participants)*8
    assert sum(sample.metadata['selected_genotype_class_counts_H_R_A']) == len(d.rows)
    return sample.metadata


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Compute allocation required; do not run on login node.')
    off = generate('nuisance_off')
    on = generate('nuisance_on')
    report = {'regimes':[verify(off), verify(on)]}
    assert off.design.participants == on.design.participants
    np.testing.assert_array_equal(off.preselection_stage_probabilities,
                                  on.preselection_stage_probabilities)
    assert [(r.individual,r.gene,r.n,r.orientation) for r in off.candidate_tags] == [
            (r.individual,r.gene,r.n,r.orientation) for r in on.candidate_tags]
    assert off.metadata['candidate_genotype_class_counts_H_R_A'] == on.metadata['candidate_genotype_class_counts_H_R_A']
    repeated = generate('nuisance_off')
    assert off.candidate_tags == repeated.candidate_tags
    assert off.selected_keys == repeated.selected_keys
    out = Path(os.environ['CODEX_REC_OUTPUT']) / ('synthetic_p1_generator_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    report.update(status='synthetic_generator_invariants_passed', fit_performed=False,
                  P1_size_or_power_evaluated=False, python=sys.version,
                  numpy=np.__version__, scipy=scipy.__version__,
                  hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in ('synthetic_p1_preascertainment.py',
                                       'check_synthetic_p1_preascertainment.py',
                                       'likelihood_v1.py', 'design_v1.py')})
    (out/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],'selected_rows':[len(off.design.rows),len(on.design.rows)]}),flush=True)


if __name__ == '__main__':
    main()
