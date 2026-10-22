"""Summarize AN5 refit diagnostics into one per-row table (login-node scale: small TSVs)."""
import pandas as pd, glob, os, numpy as np, json, sys
A = os.path.dirname(os.path.abspath(__file__))
run = sorted(d for d in glob.glob(A + '/run_2*') if os.path.isdir(d))[-1]
D = run + '/diag'
r2 = pd.read_csv(A + '/../pph4_1_rows_positive_t12.tsv', sep='\t')
plan = pd.read_csv(run + '/plan_rows.tsv', sep='\t')
out = []
for _, p in plan.iterrows():
    tag = f"{p.gene}__{p.gwas_name}"
    f = f"{D}/{tag}__fits.tsv"
    if not os.path.exists(f):
        out.append({'gene': p.gene, 'gwas': p.gwas_name, 'role': p.role, 'status': 'no diagnostic output'}); continue
    d = pd.read_csv(f, sep='\t').set_index('fit')
    k = pd.read_csv(f"{D}/{tag}__kriging_all.tsv.gz", sep='\t')
    k = k[~k.snp.duplicated()].set_index('snp')
    a = d.loc['a_primary_diag1e-3']
    hit = a.best_hit1
    zhit = k.z.get(hit, np.nan) if isinstance(hit, str) else np.nan
    zsd_hit = k.z_std_diff.get(hit, np.nan) if isinstance(hit, str) else np.nan
    # the GWAS credible set that carries the best pair
    leads = str(a.cs_lead).split(';'); sizes = str(a.cs_sizes).split(';'); lz = str(a.cs_lead_z).split(';')
    rec = {'gene': p.gene, 'gwas': p.gwas_name, 'role': p.role, 'n_snps': int(a.n_snps),
           'lambda_s': a.lambda_s, 'lead_snp_max_absz': a.lead_snp, 'lead_z': round(a.lead_z, 2),
           'kriging_flag_logLR2_absz2': int(a.n_flag_logLR2_absz2), 'max_abs_z_std_diff': round(a.max_abs_z_std_diff, 1),
           'n_abs_z_std_diff_gt3': int(a.n_abs_z_std_diff_gt3),
           'primary_converged': bool(a.converged), 'primary_niter': int(a.niter), 'primary_n_cs': int(a.n_cs),
           'primary_cs_sizes': a.cs_sizes, 'primary_cs_lead_z': a.cs_lead_z,
           'best_pair_gwas_hit1': hit, 'z_at_best_hit1': round(zhit, 2) if zhit == zhit else None,
           'z_std_diff_at_best_hit1': round(zsd_hit, 2) if zsd_hit == zsd_hit else None,
           'pp4_primary_refit': a.best_pp4}
    m = r2[(r2.gene == p.gene) & (r2.gwas_name == p.gwas_name)]
    rec['pp4_r2_release'] = float(m['PP.H4.susie'].iloc[0]) if len(m) else None
    for lab, col in [('b_diag1e-6', 'pp4_diag1e-6'), ('c_diag1e-3_estRV', 'pp4_estRV'), ('d_diag1e-3_dropKrigingFlag', 'pp4_drop_flagged')]:
        if lab in d.index and d.loc[lab].fit_ok == True:
            rec[col] = d.loc[lab].best_pp4
            rec[col.replace('pp4', 'ncs')] = d.loc[lab].n_cs
            rec[col.replace('pp4', 'conv')] = d.loc[lab].converged
        else:
            rec[col] = None
    out.append(rec)
o = pd.DataFrame(out)
o.to_csv(A + '/pph4_1_diagnostic_summary.tsv', sep='\t', index=False)
pd.set_option('display.width', 300); pd.set_option('display.max_columns', 40)
print(o.to_string())
