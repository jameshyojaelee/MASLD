#!/usr/bin/env python3
"""Require all run reads once, one capture library, and unique integer counts."""
import argparse
import json
from pathlib import Path
import os

import h5py
import numpy as np
import pandas as pd

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')


def check(directory):
    assert os.environ.get('SLURM_JOB_ID')
    spec=json.loads((directory/'input_specification.json').read_text())
    reads=[]; mapping=[]
    for run in spec['runs']:
        row=pd.read_csv(ROOT/f'data/GSE244832/cellranger/{run}/outs/metrics_summary.csv',dtype=str).iloc[0]
        reads.append(int(row['Number of Reads'].replace(',','')))
        mapping.append(float(row['Reads Mapped Confidently to Genome'].strip('%'))/100)
    row=pd.read_csv(directory/'count/outs/metrics_summary.csv',dtype=str).iloc[0]
    observed=int(row['Number of Reads'].replace(',',''))
    assert observed==sum(reads),f'Joint count read total {observed} differs from per-run total {sum(reads)}'
    weighted_mapping=float(np.average(mapping,weights=reads))
    joint_mapping=float(row['Reads Mapped Confidently to Genome'].strip('%'))/100
    assert joint_mapping >= weighted_mapping-.05,'Joint alignment lost more than five percentage points versus its input runs'
    with h5py.File(directory/'count/outs/molecule_info.h5','r') as h:
        info=h['library_info'][()]
        if isinstance(info,np.ndarray): info=b''.join(info.tolist())
        libraries=json.loads(info.decode() if isinstance(info,bytes) else info)
        assert len(libraries)==1 and libraries[0]['gem_group']==1
    with h5py.File(directory/'count/outs/filtered_feature_bc_matrix.h5','r') as h:
        m=h['matrix']; bc=m['barcodes'].asstr()[:]; data=m['data'][:]
        assert len(bc)==len(set(bc))
        assert (data>=0).all() and np.isfinite(data).all() and np.equal(data,np.floor(data)).all()
    result=dict(status='passed',library=spec['library'],runs=len(spec['runs']),joint_reads=observed,
                sum_per_run_reads=sum(reads),single_capture_library=True,unique_barcodes=True,integer_counts=True,
                joint_confident_genome_mapping=joint_mapping,weighted_run_confident_genome_mapping=weighted_mapping)
    (directory/'scientific_checks.json').write_text(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    print(json.dumps(check(p.parse_args().directory)),flush=True)
