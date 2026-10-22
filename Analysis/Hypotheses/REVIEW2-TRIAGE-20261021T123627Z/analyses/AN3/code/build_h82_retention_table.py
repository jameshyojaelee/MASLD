"""AN3 task 4: H82 joint-adjustment retention for the 11 primary-test contrasts (descriptive; PLAN D4).
Reads only saved H82 outputs; no refit. Retained = released calls that stay called (same rule) under the arm."""
import pandas as pd
R = 'Analysis/Hypotheses/H82-20261015T125052Z/run-20261015T131530Z'
OUT = 'Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/h82'
arms = pd.read_csv(f'{R}/d1_arm_table.tsv', sep='\t')
pt = pd.read_csv(f'{R}/primary_tests.tsv', sep='\t')
assert len(pt) == 11
rows = []
for _, p in pt.iterrows():
    a = arms[arms.contrast == p.contrast].set_index('arm')
    assert a.loc['base', 'calls'] == p.released == a.loc['joint', 'released']
    assert a.loc['joint', 'calls'] == p.joint_calls and a.loc['joint', 'retained'] == p.joint_retained
    r = dict(contrast=p.contrast, family=p.family, n_units=int(a.loc['base', 'n_units']), released=int(p.released))
    for arm in ['lib', 'gini', 'nf', 'joint']:
        r[f'{arm}_calls'] = int(a.loc[arm, 'calls']); r[f'{arm}_retained'] = int(a.loc[arm, 'retained'])
        r[f'{arm}_retained_frac'] = round(float(a.loc[arm, 'retained_frac']), 3)
    r.update(joint_null_median=round(p.null_median, 3), joint_null_q025=round(p.null_q025, 3),
             joint_null_q975=round(p.null_q975, 3), class_as_written=p['class'])
    rows.append(r)
t = pd.DataFrame(rows)
t.to_csv(f'{OUT}/h82_joint_adjustment_retention_11_contrasts.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(t[['contrast', 'family', 'n_units', 'released', 'lib_retained_frac', 'gini_retained_frac', 'nf_retained_frac',
         'joint_calls', 'joint_retained', 'joint_retained_frac', 'joint_null_median', 'joint_null_q025', 'joint_null_q975', 'class_as_written']].to_string(index=False))
