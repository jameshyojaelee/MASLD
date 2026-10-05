#!/usr/bin/env python3
"""Collect synthetic start receipts; boundary phases remain explicitly deferred.

This is not a complete observed-fit result. A numerical continuation must execute
the inherited conditional boundary restart/null snap and final diagnostics.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def require(ok, message):
    if not ok: raise ValueError(message)


def read(path, digest):
    require(path.stat().st_size<=32*1024**2, 'Input packet exceeds32MiB')
    body = path.read_bytes()
    require(hashlib.sha256(body).hexdigest()==digest, 'Receipt pin differs: '+str(path))
    return json.loads(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-self-sha',required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--inspect',action='store_true')
    for start in ('warm','cold'):
        parser.add_argument('--'+start+'-dir',type=Path,required=True)
        parser.add_argument('--'+start+'-protocol-sha',required=True)
        parser.add_argument('--'+start+'-terminal-sha',required=True)
        parser.add_argument('--'+start+'-terminal',choices=('start_result','failure'),required=True)
    args = parser.parse_args()
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==args.expected_self_sha,
            'Collector source pin differs')
    receipts, protocols = {}, {}
    for start in ('warm','cold'):
        directory = getattr(args,start+'_dir')
        terminal = getattr(args,start+'_terminal')
        protocols[start] = read(directory/'start_protocol.json',getattr(args,start+'_protocol_sha'))
        receipts[start] = read(directory/(terminal+'.json'),getattr(args,start+'_terminal_sha'))
        require(protocols[start]['start']==start and receipts[start]['start']==start,
                'Start label differs')
        require(protocols[start]['kappa_attempt_order']==[0.,-.1,.1]
                and protocols[start]['fixed_v'] is None,'Inherited attempt specification differs')
    for field in ('shared_null_sha','shared_protocol_sha','base_sha','model_sha','batch_size',
                  'observed_count_vector_sha256','pins','environment'):
        require(protocols['warm'][field]==protocols['cold'][field],'Shard family differs: '+field)
    failed = [s for s in ('warm','cold') if getattr(args,s+'_terminal')=='failure']
    selected = None
    if not failed:
        for start in ('warm','cold'):
            fit = receipts[start]['fit']
            require(len(fit['parameters'])==28 and len(fit['gradient'])==29
                    and all(isinstance(x,(int,float)) and math.isfinite(x)
                            for x in [fit['loglik'],fit['v'],*fit['parameters'],*fit['gradient']])
                    and fit['v']>=0 and isinstance(fit['converged'],bool),
                    'Successful start returned invalid/nonfinite fit')
        selected = ('cold' if receipts['cold']['fit']['loglik']-receipts['warm']['fit']['loglik']
                    >=1e-6 else 'warm')
    result = dict(receipts=receipts,protocols=protocols,failed_starts=failed,
        native_outer_fit_would_raise=bool(failed),selected_before_boundary_phases=selected,
        selected_fit=None if selected is None else receipts[selected]['fit'],
        start_logliks=None if failed else {s:receipts[s]['fit']['loglik'] for s in ('warm','cold')},
        start_converged=None if failed else {s:receipts[s]['fit']['converged'] for s in ('warm','cold')},
        selection_rule='cold iff coldLL-warmLL>=1e-6; convergence does not filter selection',
        whole_fit_complete=False,fit_comparison_complete=False,
        boundary_restart_deferred=True,null_snap_deferred=True,final_diagnostics_deferred=True,
        continuation_rules=['If selected v==0 and right score>1e-6, inherited single-start restart at tau_max^2; replace only for strictly greater LL.',
            'Then if selected v==0 and nullLL>=selectedLL-1e-6, recompute exact marginal at shared null; convergence=null.converged and right score<=1e-6.'],
        no_fits_calibration_or_model_validation_claimed=True)
    if args.inspect:
        print(json.dumps(dict(receipt_pins_checked=True,failed_starts=failed,
                             whole_fit_complete=False,result_written=False)))
        return
    require(args.output is not None,'Explicit exclusive output required')
    body = (json.dumps(result,allow_nan=False,indent=2)+'\n').encode()
    require(len(body)<=32*1024**2,'Collection exceeds32MiB')
    with args.output.open('xb') as handle: handle.write(body)


if __name__=='__main__': main()
