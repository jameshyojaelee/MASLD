#!/usr/bin/env python3
"""Helpers for fixing the derived copies of the Pazoki 2022 PDFF GWAS (GCST90267352, PMID 36402844).

The effect/other allele labels of the study were exchanged on disk (see
GWAS/orientation_checks/pazoki_pdff-20261010T022251Z). The MR_Data copy, which was corrected there, is the
reference here. Subcommands:

  selftest
  hg19-columns   IN_GZ OUT_GZ REPORT         set hg19_chr and hg19_pos from chromosome and base_pair_location
  positions      B E_SWAPPED E_OUT REPORT [OFFSETS_REF]
                                             restore the position column of Pazoki_PDFF_reformatted_hg19.tsv from
                                             the MR_Data hg19 copy B, row by row (E_SWAPPED already has the label
                                             exchange applied)
  lead-check     BEFORE AFTER C REPORT       check the 11-row lead table against the corrected preprocessed copy
  agreement      REPORT A B_OLD B_NEW E_OLD E_NEW R1 R2
                                             share of variants whose chromosome:position appears in R1 or R2
  update-sources IN OUT STUDY OLD_SHA NEW_SHA

No randomness is used.
"""
import gzip
import io
import json
import sys

import numpy as np
import pandas as pd

REL_TOL = 1e-9
CHROMS = [str(i) for i in range(1, 23)]
WINDOW = 5_000_000


def close(x, y):
    return abs(x - y) <= REL_TOL * max(abs(x), abs(y))


def match_subsequence(e, a, window=WINDOW):
    """For each row of E return the index of the matching row of A.

    E and A are tuples of parallel lists (chromosome, allele1, allele2, beta, se, p). E must be a subsequence
    of A in the same order, so A is scanned forward once. A row with no match within `window` rows raises.
    """
    ec, e1, e2, eb, es, ep = e
    ac, a1, a2, ab, as_, ap = a
    n, m = len(ec), len(ac)
    idx = np.empty(n, dtype=np.int64)
    j = 0
    for i in range(n):
        k = j
        stop = min(m, j + window)
        while k < stop:
            if (ac[k] == ec[i] and a1[k] == e1[i] and a2[k] == e2[i]
                    and close(ab[k], eb[i]) and close(as_[k], es[i]) and close(ap[k], ep[i])):
                break
            k += 1
        else:
            raise ValueError(f"row {i} of E has no match in the {window} rows of A after row {j}")
        idx[i] = k
        j = k + 1
    return idx


def selftest():
    a = (['1', '1', '1', '2', '2', '3'],
         ['A', 'A', 'C', 'G', 'G', 'T'],
         ['G', 'G', 'T', 'A', 'A', 'C'],
         [0.1, 0.2, 0.3, 0.4, 0.4, 0.5],
         [0.01, 0.01, 0.01, 0.01, 0.01, 0.02],
         [0.5, 0.01, 0.2, 1e-300, 1e-300, 0.3])
    # rows 0, 2, 3, 4, 5 of A, with last-digit noise on the floats
    pick = [0, 2, 3, 4, 5]
    e = tuple([col[i] for i in pick] for col in a)
    e = (e[0], e[1], e[2], [x * (1 + 1e-14) for x in e[3]], [x * (1 - 1e-14) for x in e[4]], e[5])
    got = match_subsequence(e, a, window=10)
    assert got.tolist() == pick, got.tolist()
    bad = (['1'], ['A'], ['G'], [0.999], [0.01], [0.5])
    try:
        match_subsequence(bad, a, window=10)
    except ValueError:
        pass
    else:
        raise AssertionError("a row absent from A was matched")
    # E out of order cannot be matched by a forward scan
    reordered = tuple([col[i] for i in (2, 0)] for col in a)
    try:
        match_subsequence(reordered, a, window=10)
    except ValueError:
        pass
    else:
        raise AssertionError("an out-of-order row was matched")
    print("selftest ok")


def write_report(path, obj):
    with open(path, 'w') as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)
        fh.write('\n')


def fix_hg19_columns(src, dst, report):
    """hg19_chr := chromosome, hg19_pos := base_pair_location, every other byte unchanged.

    The source positions of this study are GRCh37; the old hg19 columns were lifted a second time.
    """
    n = n_chr = n_pos = 0
    with gzip.open(src, 'rt', encoding='latin-1', newline='') as fin, open(dst, 'wb') as raw:
        gz = gzip.GzipFile(filename='', mode='wb', fileobj=raw, compresslevel=6, mtime=0)
        fout = io.TextIOWrapper(gz, encoding='latin-1', newline='\n')
        header = fin.readline().rstrip('\n').split('\t')
        ci = {h: header.index(h) for h in ('chromosome', 'base_pair_location', 'hg19_chr', 'hg19_pos')}
        assert header[-2:] == ['hg19_chr', 'hg19_pos'], header
        fout.write('\t'.join(header) + '\n')
        for line in fin:
            f = line.rstrip('\n').split('\t')
            if len(f) != len(header):
                raise ValueError(f"field count {len(f)} != {len(header)} at data row {n + 1}")
            n += 1
            if f[ci['hg19_chr']] != f[ci['chromosome']]:
                n_chr += 1
            if f[ci['hg19_pos']] != f[ci['base_pair_location']]:
                n_pos += 1
            f[ci['hg19_chr']] = f[ci['chromosome']]
            f[ci['hg19_pos']] = f[ci['base_pair_location']]
            fout.write('\t'.join(f) + '\n')
        fout.flush()
        fout.detach()
        gz.close()
    write_report(report, {'rows': n, 'hg19_chr_changed': n_chr, 'hg19_pos_changed': n_pos})
    print(f"hg19-columns: rows {n}, hg19_chr changed {n_chr}, hg19_pos changed {n_pos}")


def read_study_a(path):
    """Rows of the corrected MR_Data file on chromosomes 1 to 22, in file order."""
    df = pd.read_csv(
        path, sep='\t', usecols=['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta',
                                 'standard_error', 'p_value'],
        dtype={'chromosome': str, 'base_pair_location': np.int64, 'effect_allele': str, 'other_allele': str,
               'beta': np.float64, 'standard_error': np.float64, 'p_value': np.float64},
        keep_default_na=False, na_values=['NA'])
    return df[df['chromosome'].isin(CHROMS)].reset_index(drop=True)


def restore_positions(b_path, e_in, e_out, report, offsets_ref=None):
    """Set the position column of the reformatted copy to the source GRCh37 position.

    The reformatted copy (E) and the MR_Data hg19 copy (B) hold the same variants in the same order with the same
    second-lift positions, and B keeps the source position in base_pair_location. Each E row is therefore paired
    with the B row at the same index, and the pairing is checked on chromosome, both alleles, beta, se, p and the
    old lifted position. Any mismatch stops the run before a file is written. Row-by-row pairing is used because
    the source file repeats chromosome, alleles, beta, se and p on many rows (a search on those columns alone can
    pick the wrong twin).
    """
    b = pd.read_csv(
        b_path, sep='\t', usecols=['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta',
                                   'standard_error', 'p_value', 'hg19_pos'],
        dtype={'chromosome': str, 'base_pair_location': np.int64, 'effect_allele': str, 'other_allele': str,
               'beta': np.float64, 'standard_error': np.float64, 'p_value': np.float64, 'hg19_pos': np.int64},
        keep_default_na=False, na_values=['NA'])
    e = pd.read_csv(
        e_in, sep='\t', dtype={'chromosome': str, 'position': np.int64, 'allele1': str, 'allele2': str,
                               'beta': np.float64, 'se': np.float64, 'pval': np.float64},
        keep_default_na=False, na_values=['NA'])
    assert list(e.columns) == ['chromosome', 'position', 'allele1', 'allele2', 'beta', 'se', 'pval'], list(e.columns)
    print(f"B rows: {len(b)}; E rows: {len(e)}", flush=True)
    out = {'rows_B': len(b), 'rows_E': len(e)}
    if len(b) != len(e):
        write_report(report, out)
        raise SystemExit("B and E have different row counts")

    def same(x, y):
        return np.isclose(x.to_numpy(), y.to_numpy(), rtol=REL_TOL, atol=0.0)

    agree = {'chromosome': (e['chromosome'] == b['chromosome']).to_numpy(),
             'allele1_vs_effect_allele': (e['allele1'] == b['effect_allele']).to_numpy(),
             'allele2_vs_other_allele': (e['allele2'] == b['other_allele']).to_numpy(),
             'beta': same(e['beta'], b['beta']), 'se': same(e['se'], b['standard_error']),
             'pval': same(e['pval'], b['p_value']),
             'old_position_vs_old_hg19_pos': (e['position'] == b['hg19_pos']).to_numpy()}
    out['pairing_mismatches'] = {k: int((~v).sum()) for k, v in agree.items()}
    keys = ['chromosome', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value']
    dup = b.duplicated(subset=keys, keep=False)
    out['rows_sharing_chromosome_alleles_beta_se_p_with_another_row'] = int(dup.sum())
    out['of_which_also_sharing_position'] = int(b[dup].duplicated(subset=keys + ['base_pair_location'],
                                                                  keep=False).sum())
    if any(out['pairing_mismatches'].values()):
        write_report(report, out)
        print(json.dumps(out, indent=1))
        raise SystemExit("row pairing between E and B failed")
    new_pos = b['base_pair_location'].to_numpy()
    offset = e['position'].to_numpy() - new_pos
    changed = int((offset != 0).sum())
    medians = pd.Series(offset).groupby(e['chromosome'].to_numpy()).median()
    out.update({'position_changed': changed, 'position_unchanged': int(len(e) - changed),
                'median_offset_old_minus_new_by_chromosome': {str(k): float(v) for k, v in medians.items()}})
    if offsets_ref:
        ref = pd.read_csv(offsets_ref, sep='\t', dtype={'chromosome': str}).set_index('chromosome')['median']
        out['max_abs_median_diff_vs_MR_hg19_offsets'] = float((medians - ref.reindex(medians.index)).abs().max())
    with open(e_in, encoding='latin-1') as fin, open(e_out, 'w', encoding='latin-1', newline='\n') as fout:
        header = fin.readline()
        fout.write(header)
        pos_col = header.rstrip('\n').split('\t').index('position')
        for i, line in enumerate(fin):
            f = line.rstrip('\n').split('\t')
            f[pos_col] = str(new_pos[i])
            fout.write('\t'.join(f) + '\n')
    assert i + 1 == len(e), (i + 1, len(e))
    write_report(report, out)
    print(json.dumps({k: v for k, v in out.items() if k != 'median_offset_old_minus_new_by_chromosome'}))


def lead_check(before, after, c_path, report):
    lb = pd.read_csv(before, sep='\t', dtype=str, keep_default_na=False)
    la = pd.read_csv(after, sep='\t', dtype=str, keep_default_na=False)
    cols = ['chromosome', 'position', 'locus', 'allele1', 'allele2', 'beta', 'se', 'pval']
    assert list(lb.columns) == cols and list(la.columns) == cols, (list(lb.columns), list(la.columns))
    swapped = bool((lb['allele1'] == la['allele2']).all() and (lb['allele2'] == la['allele1']).all()
                   and (lb[['chromosome', 'position', 'locus', 'beta', 'se', 'pval']]
                        == la[['chromosome', 'position', 'locus', 'beta', 'se', 'pval']]).all().all())
    c = pd.read_csv(c_path, sep='\t', usecols=['chromosome', 'position', 'allele1', 'allele2', 'beta'],
                    dtype={'chromosome': str, 'position': str, 'allele1': str, 'allele2': str, 'beta': np.float64})
    keys = set(zip(la['chromosome'], la['position']))
    c = c[[k in keys for k in zip(c['chromosome'], c['position'])]]
    rows = []
    for _, r in la.iterrows():
        hit = c[(c['chromosome'] == r['chromosome']) & (c['position'] == r['position'])
                & (c['allele1'] == r['allele1']) & (c['allele2'] == r['allele2'])]
        ok = len(hit) == 1 and close(float(hit['beta'].iloc[0]), float(r['beta']))
        rows.append({'chromosome': r['chromosome'], 'position': r['position'], 'locus': r['locus'],
                     'allele1': r['allele1'], 'allele2': r['allele2'], 'matching_rows_in_corrected_C': len(hit),
                     'beta_equal': bool(ok)})
    res = pd.DataFrame(rows)
    res.to_csv(report, sep='\t', index=False)
    print(res.to_string(index=False))
    if not swapped or not res['beta_equal'].all():
        raise SystemExit("lead-check failed")
    print(f"lead-check ok: {len(res)} rows, labels exchanged, every row equals the corrected C row")


def position_keys(chrom, pos):
    chrom = pd.Series(chrom).astype(str)
    keep = chrom.isin(CHROMS).to_numpy()
    return chrom[keep].astype(np.int64).to_numpy() * 10**9 + np.asarray(pos)[keep].astype(np.int64)


def agreement(report, a_path, b_old, b_new, e_old, e_new, r1, r2):
    ref = []
    for p in (r1, r2):
        d = pd.read_csv(p, sep='\t', usecols=['chromosome', 'position'], dtype={'chromosome': str, 'position': np.int64})
        ref.append(position_keys(d['chromosome'], d['position']))
    ref = np.unique(np.concatenate(ref))
    out = {'reference_variants_R1_R2': int(len(ref))}
    a = read_study_a(a_path)
    sets = {'A_base_pair_location': (a['chromosome'], a['base_pair_location'])}
    for name, p in (('B_old_hg19', b_old), ('B_new_hg19', b_new)):
        d = pd.read_csv(p, sep='\t', usecols=['hg19_chr', 'hg19_pos'], dtype={'hg19_chr': str, 'hg19_pos': np.int64})
        sets[name] = (d['hg19_chr'], d['hg19_pos'])
    for name, p in (('E_old', e_old), ('E_new', e_new)):
        d = pd.read_csv(p, sep='\t', usecols=['chromosome', 'position'], dtype={'chromosome': str, 'position': np.int64})
        sets[name] = (d['chromosome'], d['position'])
    for name, (c, pos) in sets.items():
        k = position_keys(c, pos)
        out[name] = {'variants': int(len(k)), 'share_in_R1_or_R2': float(np.isin(k, ref).mean())}
    write_report(report, out)
    print(json.dumps(out, indent=1))
    share = {k: v['share_in_R1_or_R2'] for k, v in out.items() if isinstance(v, dict)}
    a_share = share['A_base_pair_location']
    ok = (a_share > 0.5 and abs(share["E_new"] - a_share) < 0.02 and abs(share["B_new_hg19"] - a_share) < 0.02
          and share['E_old'] < 0.1 and share['B_old_hg19'] < 0.1)
    if not ok:
        raise SystemExit("position agreement check failed")
    print("agreement ok")


def update_sources(src, dst, study, old_sha, new_sha):
    with open(src, encoding='latin-1') as fh:
        lines = fh.read().split('\n')
    header = lines[0].split('\t')
    i_name, i_sha = header.index('study_name'), header.index('sha256')
    n = 0
    for k in range(1, len(lines)):
        f = lines[k].split('\t')
        if len(f) > i_sha and f[i_name] == study:
            assert f[i_sha] == old_sha, (f[i_sha], old_sha)
            f[i_sha] = new_sha
            lines[k] = '\t'.join(f)
            n += 1
    assert n == 1, n
    with open(dst, 'w', encoding='latin-1', newline='\n') as fh:
        fh.write('\n'.join(lines))
    print(f"update-sources: {study} sha256 {old_sha[:8]} -> {new_sha[:8]}")


def main(argv):
    cmd, args = argv[1], argv[2:]
    if cmd == 'selftest':
        selftest()
    elif cmd == 'hg19-columns':
        fix_hg19_columns(*args)
    elif cmd == 'positions':
        restore_positions(*args)
    elif cmd == 'lead-check':
        lead_check(*args)
    elif cmd == 'agreement':
        agreement(*args)
    elif cmd == 'update-sources':
        update_sources(*args)
    else:
        raise SystemExit(__doc__)


if __name__ == '__main__':
    main(sys.argv)
