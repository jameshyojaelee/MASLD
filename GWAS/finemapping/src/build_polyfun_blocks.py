#!/usr/bin/env python
"""
build_polyfun_blocks.py — convert PolyFun precomputed UKBB LD (.npz + .gz pairs,
3 Mb windows tiled every 1 Mb) into per-Berisa-Pickrell-block .bim + .ld files
matching the existing finemapping pipeline convention.

Phase 9 rewrite (2026-04-28):
  - Previous version assumed PolyFun windows are non-overlapping 3 Mb tiles
    and `bdiag()`-glued multiple windows when a BP block straddled a 3 Mb
    edge. Cross-window LD entries went to zero, which broke SuSiE-CS:
    `coloc.susie` returned NULL for every eGene, all 374 SuSiE-COLOC tasks
    silently fell back to ABF.
  - Reality: PolyFun windows are 3 Mb wide placed every 1 Mb (overlap 2 Mb).
    96.6% of BP blocks (≤ 3 Mb wide) fit fully inside at least one window;
    we pick the smallest such window and use its full intra-block LD. The
    remaining 3.4% of blocks (centromere/pericentromeric, > 3 Mb wide) get
    truncated to the window with maximum overlap.
  - Result: per-block LD is now full-rank UKBB-quality, structurally
    equivalent to sghatan v1 / 1KG (PLINK `--r square`) — no `bdiag`.

Input:
  data/ld_ref/polyfun_eur/raw/chr<N>_<start>_<end>.{npz,gz}  (2.85 TB)
    — start advances every 1_000_000; window width is 3_000_000.
  data/ld_ref/1kg_eur/approx_LD_blocks.txt  (1703 BP blocks)

Output:
  data/ld_ref/polyfun_eur/chr<N>/<bs>.<be>/<bs>.<be>.bim
  data/ld_ref/polyfun_eur/chr<N>/<bs>.<be>/<bs>.<be>.ld

Usage:
  Single chr:    python build_polyfun_blocks.py 22
  Sbatch array:  SLURM_ARRAY_TASK_ID dispatches chr 1..22 (one chr per task).

Format details (PolyFun finemapper.py wiki):
  .npz: scipy.sparse lower-triangular correlation matrix; load + symmetrize.
  .gz : whitespace-separated rsid / chromosome / position / allele1 / allele2.
  Variant order in .gz matches row/col order in .npz.
"""
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

BASE = Path('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping')
RAW_DIR = BASE / 'data/ld_ref/polyfun_eur/raw'
OUT_DIR = BASE / 'data/ld_ref/polyfun_eur'
BLOCKS_FILE = BASE / 'data/ld_ref/1kg_eur/approx_LD_blocks.txt'
WIN_WIDTH = 3_000_000
WIN_STRIDE = 1_000_000


def chr_from_env_or_arg():
    """Get chr number from SLURM_ARRAY_TASK_ID or argv[1]."""
    raw = os.environ.get('SLURM_ARRAY_TASK_ID')
    if raw is None and len(sys.argv) > 1:
        raw = sys.argv[1]
    if raw is None:
        raise RuntimeError('Need SLURM_ARRAY_TASK_ID or argv[1] for chr number')
    return int(raw)


def load_polyfun_window(npz_path, gz_path):
    """Load a 3 Mb PolyFun window: returns (R, variants_df).
    R is the symmetric NxN correlation matrix; variants_df is in row order.
    """
    R_lower = sparse.load_npz(str(npz_path)).toarray().astype(np.float32)
    R = R_lower + R_lower.T
    np.fill_diagonal(R, 1.0)
    assert R.shape[0] == R.shape[1], f'{npz_path} not square'

    df = pd.read_csv(gz_path, sep=r'\s+', compression='gzip', engine='python')
    df.columns = [c.lower() for c in df.columns]
    needed = ['rsid', 'chromosome', 'position', 'allele1', 'allele2']
    for c in needed:
        if c not in df.columns:
            raise RuntimeError(f'{gz_path} missing column {c}; got {list(df.columns)}')
    df = df[needed]

    if len(df) != R.shape[0]:
        raise RuntimeError(
            f'mismatch {gz_path}: {len(df)} variants vs {R.shape[0]} R rows'
        )
    return R, df


def list_existing_window_starts(chrom):
    """Return sorted list of window start positions on disk for a chromosome."""
    starts = []
    for p in RAW_DIR.glob(f'chr{chrom}_*.npz'):
        try:
            # filename: chr<N>_<start>_<end>.npz
            parts = p.stem.split('_')
            start = int(parts[1])
            starts.append(start)
        except Exception:
            continue
    return sorted(set(starts))


def find_best_window_for_block(chrom, block_start, block_end, win_starts):
    """Return (win_start, win_end, status, eff_start, eff_end) — see docstring.

    Strategy:
      1. If block_width <= WIN_WIDTH and a window FULLY contains [block_start,
         block_end), pick the window with smallest start (deterministic). Use
         the block's exact range as the effective subset. Status = 'FIT'.
      2. Otherwise (block > WIN_WIDTH, or no fully-containing window exists at
         a chromosome edge), pick the window with the largest overlap to the
         block. Truncate effective range to the intersection. Status =
         'TRUNCATED'.
      3. If no window overlaps the block at all, return None. Status =
         'NO_WINDOW'.
    """
    # Filter to windows that touch the block (overlap > 0)
    candidates = []
    for ws in win_starts:
        we = ws + WIN_WIDTH
        if we <= block_start or ws >= block_end:
            continue
        candidates.append((ws, we))
    if not candidates:
        return (None, None, 'NO_WINDOW', None, None)

    block_width = block_end - block_start
    if block_width <= WIN_WIDTH:
        # Containing windows: ws <= block_start AND we >= block_end
        contains = [(ws, we) for ws, we in candidates
                    if ws <= block_start and we >= block_end]
        if contains:
            ws, we = min(contains)  # smallest start (deterministic)
            return (ws, we, 'FIT', block_start, block_end)
        # No fully-containing window: fall through to overlap-max truncation

    # Overlap-max truncation (wide blocks or chromosome edges)
    best = None
    for ws, we in candidates:
        ovl = min(we, block_end) - max(ws, block_start)
        if best is None or ovl > best[2]:
            best = (ws, we, ovl)
    ws, we, _ = best
    eff_start = max(ws, block_start)
    eff_end = min(we, block_end)
    return (ws, we, 'TRUNCATED', eff_start, eff_end)


def build_block(chrom, block_start, block_end, win_starts, win_cache, stats):
    """Build .bim + .ld for a single block. Returns (n_variants, status_str)."""
    out_subdir = OUT_DIR / f'chr{chrom}' / f'{block_start}.{block_end}'
    bim_path = out_subdir / f'{block_start}.{block_end}.bim'
    ld_path = out_subdir / f'{block_start}.{block_end}.ld'

    # Idempotent: skip if both outputs exist + non-empty + bim row count
    # matches expected ld dimension.
    if bim_path.exists() and ld_path.exists() and bim_path.stat().st_size > 0:
        try:
            n_bim = sum(1 for _ in open(bim_path))
            with open(ld_path) as fh:
                first = fh.readline()
                n_ld_cols = len(first.split())
            if n_bim == n_ld_cols and n_bim > 0:
                return n_bim, 'SKIP_EXISTS'
        except Exception:
            pass  # rebuild

    win_start, win_end, status, eff_start, eff_end = find_best_window_for_block(
        chrom, block_start, block_end, win_starts)

    if status == 'NO_WINDOW':
        return 0, 'NO_WINDOW'

    cache_key = f'{win_start}_{win_end}'
    if cache_key not in win_cache:
        npz_path = RAW_DIR / f'chr{chrom}_{win_start}_{win_end}.npz'
        gz_path = RAW_DIR / f'chr{chrom}_{win_start}_{win_end}.gz'
        if not (npz_path.exists() and gz_path.exists()):
            return 0, 'NO_WINDOW'
        win_cache[cache_key] = load_polyfun_window(npz_path, gz_path)
        # Cache LRU: keep at most 3 windows (a chr typically has dozens; we
        # process blocks in order so adjacent blocks share windows).
        if len(win_cache) > 3:
            first_key = next(iter(win_cache))
            if first_key != cache_key:
                del win_cache[first_key]

    R, vdf = win_cache[cache_key]

    # Subset variants to effective block range
    pos = vdf['position'].values
    in_range = (pos >= eff_start) & (pos < eff_end)
    if not in_range.any():
        return 0, 'EMPTY_BLOCK'

    idx = np.where(in_range)[0]
    bim_df = vdf.iloc[idx].copy()
    R_sub = R[np.ix_(idx, idx)]  # full intra-window LD subset to block
    n = len(bim_df)

    out_subdir.mkdir(parents=True, exist_ok=True)

    # Write .bim: chr  rsid  0  pos  allele1  allele2
    bim_out = bim_df[['chromosome', 'rsid', 'position', 'allele1', 'allele2']].copy()
    bim_out.insert(2, 'cm', 0)
    bim_out = bim_out[['chromosome', 'rsid', 'cm', 'position', 'allele1', 'allele2']]
    bim_out.to_csv(bim_path, sep='\t', header=False, index=False)

    # Write .ld as whitespace-separated dense matrix (4-decimal precision).
    # PolyFun .npz uses float32; SuSiE tolerance is ~1e-6; %.4f is sufficient.
    np.savetxt(ld_path, R_sub, fmt='%.4f', delimiter=' ')

    if status == 'TRUNCATED':
        stats['truncated'].append(
            (block_start, block_end, eff_start, eff_end, n))

    return n, status


def main():
    chrom = chr_from_env_or_arg()
    print(f'[build_polyfun_blocks] chr{chrom} start  '
          f'{time.strftime("%Y-%m-%d %H:%M:%S")}')

    blocks = pd.read_csv(BLOCKS_FILE, sep=r'\s+', engine='python')
    blocks_chr = blocks[blocks['chr'] == chrom].copy()
    print(f'  {len(blocks_chr)} BP blocks on chr{chrom}')

    win_starts = list_existing_window_starts(chrom)
    if not win_starts:
        raise RuntimeError(f'No PolyFun windows found for chr{chrom} in {RAW_DIR}')
    print(f'  {len(win_starts)} PolyFun windows on chr{chrom} '
          f'(start range {win_starts[0]}-{win_starts[-1]})')

    win_cache = {}
    stats = {'truncated': []}
    n_fit = n_trunc = n_skip = n_empty = n_nowin = 0
    t0 = time.time()
    for _, row in blocks_chr.iterrows():
        bs = int(row['start'])
        be = int(row['stop'])
        nvar, status = build_block(chrom, bs, be, win_starts, win_cache, stats)
        if status == 'FIT':
            n_fit += 1
        elif status == 'TRUNCATED':
            n_trunc += 1
        elif status == 'SKIP_EXISTS':
            n_skip += 1
        elif status == 'EMPTY_BLOCK':
            n_empty += 1
            print(f'  [empty] block {bs}.{be}: 0 variants in chosen window range')
        elif status == 'NO_WINDOW':
            n_nowin += 1
            print(f'  [no-win] block {bs}.{be}: no overlapping PolyFun window')
        if (n_fit + n_trunc) > 0 and (n_fit + n_trunc) % 25 == 0:
            elapsed = time.time() - t0
            done = n_fit + n_trunc + n_skip + n_empty + n_nowin
            print(f'  [{done}/{len(blocks_chr)}] block {bs}.{be}: {nvar} vars  ({elapsed:.1f}s)')

    print(f'\n[build_polyfun_blocks] chr{chrom} done  ({time.time() - t0:.1f}s)')
    print(f'  FIT (full containment):     {n_fit}')
    print(f'  TRUNCATED (>3Mb or edge):   {n_trunc}')
    print(f'  SKIP_EXISTS (idempotent):   {n_skip}')
    print(f'  EMPTY (0 variants):         {n_empty}')
    print(f'  NO_WINDOW (uncovered):      {n_nowin}')
    if stats['truncated']:
        print(f'\n  Truncated block details (first 10):')
        for bs, be, es, ee, n in stats['truncated'][:10]:
            print(f'    block {bs}-{be} ({be-bs:,}bp wide)  →  '
                  f'eff {es}-{ee} ({ee-es:,}bp, {n} vars)')


if __name__ == '__main__':
    main()
