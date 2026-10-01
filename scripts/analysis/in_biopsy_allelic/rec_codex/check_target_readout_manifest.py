"""Check fixed development/pilot population identity without scoring or outcomes."""
import json
import os
from pathlib import Path

import pandas as pd

from target_readout_pilot import FIELDS, read_manifest


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Compute allocation required.')
    root=Path(os.environ['CODEX_REC_OUTPUT'])
    source=root/'target_region_manifest_21995444'
    full_path=source/'development_target_manifest_1mb.tsv.gz'
    pilot_path=source/'target_readout_pilot64.tsv'
    full,full_sha=read_manifest(full_path,True)
    pilot,pilot_sha=read_manifest(pilot_path)
    matched=full.set_index('key',drop=False).loc[pilot.key].reset_index(drop=True)
    pd.testing.assert_frame_equal(matched[FIELDS[1:]],pilot[FIELDS[1:]].reset_index(drop=True))
    # Neither mode can silently accept the other's population.
    for path,mode in [(full_path,False),(pilot_path,True)]:
        try:
            read_manifest(path,mode)
        except RuntimeError:
            pass
        else:
            raise AssertionError('Population-mode mismatch accepted.')
    out=root/('target_readout_manifest_check_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    report={'status':'fixed_readout_population_checks_passed',
            'development_rows':len(full),'pilot_rows':len(pilot),
            'manifest_sha256':full_sha,'pilot_sha256':pilot_sha,
            'fold_counts':full.heldout_fold.value_counts().sort_index().to_dict(),
            'pilot_input_identity':'Exact common metadata fields equal full-manifest rows.',
            'outcomes_read':False,'inference_run':False}
    (out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
