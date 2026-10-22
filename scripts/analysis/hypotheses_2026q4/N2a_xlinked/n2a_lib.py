"""Shared N2a rules (PRESPEC_N2a.md): heterozygous call, informative site, cross-fit contribution."""
import gzip

MIN_DEPTH, MIN_MINOR, MIN_MINOR_FRAC, MAX_OTHER_FRAC = 20, 5, 0.10, 0.05
BIN_EDGES = [20, 32, 64, 128, 256, 10 ** 9]


def load_counts(path):
    """-> {(chrom, pos): (ref_p0, alt_p0, ref_p1, alt_p1, other_p0, other_p1)}"""
    d = {}
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            d[(f[0], int(f[1]))] = tuple(int(x) for x in f[2:8])
    return d


def depth_bin(n):
    for i in range(len(BIN_EDGES) - 1):
        if BIN_EDGES[i] <= n < BIN_EDGES[i + 1]:
            return i
    return None


def site_record(a, b):
    """a = counts in biopsy 1, b = counts in biopsy 2 at one site (6-tuples as in load_counts).
    Returns None unless the site is informative: depth >= 20 in each biopsy, pooled minor allele >= 5 reads and
    >= 10% of pooled REF+ALT, other bases <= 5% of all reads.
    Otherwise returns (depth_bin, x1, n1, x2, n2), the cross-fit contribution for direction 1 to 2:
    the phase allele is the majority allele among biopsy-1 parity-0 reads; x1/n1 are its count and the depth among
    biopsy-1 parity-1 reads; x2/n2 are the same for all biopsy-2 reads. A parity-0 tie gives zeros."""
    r1, a1, r2, a2 = a[0] + a[2], a[1] + a[3], b[0] + b[2], b[1] + b[3]
    n1, n2 = r1 + a1, r2 + a2
    if n1 < MIN_DEPTH or n2 < MIN_DEPTH:
        return None
    ref, alt = r1 + r2, a1 + a2
    minor = min(ref, alt)
    other = a[4] + a[5] + b[4] + b[5]
    if minor < MIN_MINOR or minor / (ref + alt) < MIN_MINOR_FRAC or other / (ref + alt + other) > MAX_OTHER_FRAC:
        return None
    if a[0] == a[1]:
        return (depth_bin(min(n1, n2)), 0, 0, 0, 0)
    phase_ref = a[0] > a[1]
    return (depth_bin(min(n1, n2)), a[2] if phase_ref else a[3], a[2] + a[3], r2 if phase_ref else a2, n2)
