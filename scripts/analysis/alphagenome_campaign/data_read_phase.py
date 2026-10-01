#!/usr/bin/env python3
"""Targeted ATAC fragment co-coverage census, not genotype or interaction validation."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil

ROOT = Path(__file__).resolve().parents[3]
F1 = ROOT / 'GWAS/finemapping/results/alphagenome_program/f1-phase-20260914T193608Z/tables/f1_read_backed_phase_candidates.tsv'
ASE = ROOT / 'GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables'
FASTA = Path('/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa')
MAPQ, BASEQ = 30, 20
FIXED = [
    ('GSE281367', 'Z01', 'chr8:9325848:A:G', 'chr8:9326086:A:G'),
    ('GSE281367', 'Z04', 'chr8:9325848:A:G', 'chr8:9326086:A:G'),
    ('GSE281367', 'Z05', 'chr8:9325848:A:G', 'chr8:9326086:A:G'),
    ('GSE281367', 'Z08', 'chr8:9325848:A:G', 'chr8:9326086:A:G'),
    ('GSE244832', 'D15', 'chr19:45326124:G:A', 'chr19:45326536:G:C'),
    ('GSE244832', 'D17', 'chr22:43928847:C:G', 'chr22:43928850:C:T'),
]


def rows(path, sep='\t'):
    op = gzip.open if path.suffix == '.gz' else open
    with op(path, 'rt', newline='') as f:
        yield from csv.DictReader(f, delimiter=sep)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def write_json(path, x):
    with path.open('x') as f: json.dump(x, f, indent=2); f.write('\n')


def table(path, data):
    with path.open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter='\t')
        w.writeheader(); w.writerows(data)


def variant(uid):
    c, p, ref, alt = uid.split(':')
    if ref == alt or len(ref) != 1 or len(alt) != 1 or set(ref+alt) - set('ACGT'):
        raise ValueError('Fixed census requires exact distinct SNV alleles')
    return c, int(p)-1, ref, alt


def call_state(values):
    return next(iter(values)) if len(values) == 1 else ('conflict' if values else None)


def molecule_state(templates):
    """Never join single-site calls from different templates to invent co-coverage."""
    all_calls = [set(), set()]
    direct = set()
    same_read = False
    for t in templates:
        for i in (0, 1): all_calls[i].update(t['calls'][i])
        states = [call_state(x) for x in t['calls']]
        if all(x in ('0', '1') for x in states): direct.add(''.join(states))
        same_read |= t['same_read']
    states = [call_state(x) for x in all_calls]
    conflict = 'conflict' in states or len(direct) > 1
    hap = next(iter(direct)) if len(direct) == 1 and not conflict else None
    return states, hap, same_read, conflict


def read_calls(read, variants):
    calls, stats, diagnostics = {}, Counter(), []
    positions = {v[1]: i for i, v in enumerate(variants)}
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos not in positions: continue
        i = positions[rpos]
        if qpos is None:
            stats['site_'+str(i+1)+'_deletion_or_skip'] += 1; continue
        if read.query_qualities is None or read.query_qualities[qpos] < BASEQ:
            stats['site_'+str(i+1)+'_low_base_quality'] += 1; continue
        # SAM/pysam SEQ is already oriented to the aligned reference strand.
        # Reverse-complementing again would corrupt reverse-read allele calls.
        base = read.query_sequence[qpos].upper()
        allele = '0' if base == variants[i][2] else '1' if base == variants[i][3] else 'other'
        calls[i] = allele
        diagnostics.append((i+1, allele, read.query_qualities[qpos]))
    return calls, stats, diagnostics


def self_tests(pysam):
    def t(a=(), b=(), same=False): return dict(calls=[set(a), set(b)], same_read=same)
    assert molecule_state([t('0', '1', True)])[1] == '01'
    assert molecule_state([t('0', '1'), t('0', '1')])[1] == '01'
    assert molecule_state([t('0'), t(b='1')])[1] is None
    assert molecule_state([t('0', '1'), t('1', '1')])[1] is None
    assert molecule_state([t('01', '1')])[3]
    assert molecule_state([t('0', 'other')])[1] is None
    assert molecule_state([t('0', '0')])[1] == '00'
    assert molecule_state([t('1', '1')])[1] == '11'
    read = pysam.AlignedSegment()
    read.query_name = 'synthetic'
    read.reference_start = 100
    read.query_sequence = 'AG'
    read.query_qualities = [40, 40]
    read.cigartuples = [(0, 2)]
    variants = [('chr1',100,'A','T'), ('chr1',101,'C','G')]
    assert read_calls(read, variants)[0] == {0:'0', 1:'1'}
    read.flag = 16
    assert read_calls(read, variants)[0] == {0:'0', 1:'1'}
    read.query_qualities = [19, 40]
    assert read_calls(read, variants)[0] == {1:'1'}
    read.query_qualities = [40, 40]
    read.cigartuples = [(0,1),(2,1),(0,1)]
    assert read_calls(read, variants)[0] == {0:'0'}
    assert read_calls(read, variants)[1]['site_2_deletion_or_skip'] == 1
    read.cigartuples = [(0,1),(3,1),(0,1)]
    assert read_calls(read, variants)[0] == {0:'0'}
    return 'template_co_coverage_duplicate_conflict_four_haplotypes_CIGAR_deletion_skip_basequality_reverse_orientation_passed'


def source_manifest():
    old = {(r['dataset'], r['site_pair']): r for r in rows(F1)}
    source_paths = {'GSE281367': ROOT/'data/GSE281367/metadata/donor_pairing.csv',
                    'GSE244832': ROOT/'data/GSE244832/metadata/donor_pairing.csv'}
    donor_maps = {s: {r['donor_id']: r for r in rows(p, ',')} for s, p in source_paths.items()}
    identities = {}
    for source, name in [('GSE281367', 'allelic_donor_counts.tsv.gz'), ('GSE244832', 'gse244832_allelic_donor_counts.tsv.gz')]:
        identities[source] = {(r['donor'], r['uid']) for r in rows(ASE/name)}
    out = []
    for source, donor, uid1, uid2 in FIXED:
        v1, v2 = variant(uid1), variant(uid2)
        if v1[0] != v2[0] or not 0 <= v1[1] < v2[1]: raise ValueError('Fixed coordinate identity differs')
        oldrow = old[(source, uid1+' | '+uid2)]
        if donor not in oldrow['donors_het_at_both'].split(','): raise ValueError('Donor absent from exact historical pair')
        if not all((donor, uid) in identities[source] for uid in [uid1, uid2]):
            raise ValueError('Exact donor/REF/ALT UID missing in source allelic table; no position-only substitution')
        d = donor_maps[source][donor]
        bam = (ROOT/f'Analysis/ATAC/Human_External/cellranger/{donor}/outs/possorted_bam.bam' if source == 'GSE281367'
               else ROOT/f'Analysis/ATAC/Human_Multiome/results/alignment/{donor}.filtered.bam')
        if not bam.is_file() or not Path(str(bam)+'.bai').is_file(): raise ValueError('Fixed donor indexed BAM missing')
        out.append(dict(source=source, donor=donor, variant1=uid1, variant2=uid2,
            source_srr=d['atac_srr'], source_gsm=d['atac_gsm'], source_atac_sample=d.get('atac_sample', donor),
            bam=str(bam), donor_map=str(source_paths[source]), exact_allelic_UID_join=True,
            source_genotype_independently_verified=False))
    return out


def scan(item, fasta, pysam):
    v = [variant(item['variant1']), variant(item['variant2'])]
    c, lo, hi = v[0][0], v[0][1], v[1][1]
    for chrom, pos, ref, alt in v:
        if fasta.fetch(chrom, pos, pos+1).upper() != ref: raise ValueError('Reference allele mismatch in fixed target')
    templates = {}
    stats, read_diagnostics = Counter(), Counter()
    bam_path = Path(item['bam'])
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        if bam.get_reference_length(c) != fasta.get_reference_length(c): raise ValueError('BAM/reference contig length mismatch')
        header = bam.header.to_dict()
        sm = {r['SM'] for r in header.get('RG', []) if r.get('SM')}
        expected = {item['donor'], item['source_srr'], item['source_gsm'], item['source_atac_sample']}
        if sm - expected:
            raise ValueError('Unexpected BAM read-group sample name; requires metadata resolution before donor counting: '+str(sm))
        for read in bam.fetch(c, lo, hi+1):
            stats['alignments_fetched'] += 1
            if stats['alignments_fetched'] > 1000000: raise ValueError('Bounded interval scan exceeds one million alignments')
            bad = []
            for field in ['is_unmapped', 'is_secondary', 'is_supplementary', 'is_qcfail', 'is_duplicate']:
                if getattr(read, field): bad.append(field)
            if not read.is_paired or not read.is_proper_pair or read.mate_is_unmapped: bad.append('not_proper_mapped_pair')
            if read.mapping_quality < MAPQ: bad.append('low_MAPQ')
            for reason in bad: stats['excluded_'+reason] += 1
            if bad: continue
            if read.next_reference_id != read.reference_id or read.template_length == 0:
                stats['excluded_invalid_template_geometry'] += 1; continue
            if item['source'] == 'GSE281367':
                barcode = read.get_tag('CB') if read.has_tag('CB') else ''
                valid = re.fullmatch(r'[ACGTN]+(?:-[0-9]+)?', barcode)
            else:
                barcode = read.query_name.split(':', 1)[0] if ':' in read.query_name else ''
                valid = re.fullmatch(r'[ACGTN]+', barcode)
            if not valid:
                stats['excluded_missing_or_unrecognized_barcode'] += 1; continue
            left = min(read.reference_start, read.next_reference_start)
            right = left + abs(read.template_length)
            if left < 0 or right <= left or read.reference_end > right:
                stats['excluded_invalid_template_geometry'] += 1; continue
            rg = read.get_tag('RG') if read.has_tag('RG') else ''
            key = (rg, read.query_name)
            geometry = (barcode, c, left, right)
            t = templates.setdefault(key, dict(geometry=geometry, calls=[set(), set()], same_read=False,
                records=0, geometry_conflict=False, read_arms=set()))
            t['records'] += 1
            arm = 1 if read.is_read1 else 2 if read.is_read2 else 0
            if arm in t['read_arms']: t['geometry_conflict'] = True
            t['read_arms'].add(arm)
            if t['geometry'] != geometry: t['geometry_conflict'] = True
            calls, call_stats, diagnostic = read_calls(read, v)
            stats.update(call_stats)
            for i, allele in calls.items(): t['calls'][i].add(allele)
            for i, allele, baseq in diagnostic:
                read_diagnostics[(i, allele, '-' if read.is_reverse else '+', arm,
                    read.mapping_quality, baseq)] += 1
            if len(calls) == 2 and all(a in ('0', '1') for a in calls.values()): t['same_read'] = True
            stats['alignments_after_filters'] += 1
    buckets = defaultdict(list)
    for t in templates.values():
        if t['geometry_conflict']:
            stats['excluded_template_identity_or_geometry_conflict'] += 1; continue
        buckets[t['geometry']].append(t)
    result, family_sizes = Counter(), Counter()
    for geometry, records in buckets.items():
        family_sizes[len(records)] += 1
        calls, hap, same_read, conflict = molecule_state(records)
        result['unique_barcode_fragment_geometries'] += 1
        if geometry[2] <= lo and geometry[3] > hi: result['geometrically_spans_both_sites_among_fetched_templates'] += 1
        any_valid = any(a in ('0', '1') for a in calls)
        if any_valid: result['molecules_with_any_callable_target_allele'] += 1
        for i, allele in enumerate(calls, 1):
            result['site_'+str(i)+'_'+str(allele)] += 1
        if conflict: result['molecules_with_call_conflict'] += 1
        if hap:
            result['haplotype_'+hap] += 1
            result['molecules_directly_cocovering_both'] += 1
            result['same_read_cocoverage' if same_read else 'mate_combined_cocoverage'] += 1
    n = result['molecules_directly_cocovering_both']
    denominator = result['molecules_with_any_callable_target_allele']
    spanning = result['geometrically_spans_both_sites_among_fetched_templates']
    for key in ['haplotype_00','haplotype_01','haplotype_10','haplotype_11']:
        result.setdefault(key, 0)
    return dict(identity=item, bam_header=header,
        BAM_sample_identity_state='header_SM_matches_curated_identifiers' if sm else 'SM_absent_donor_assignment_from_curated_path_and_source_metadata',
        bam_bytes=bam_path.stat().st_size,
        bam_mtime_ns=bam_path.stat().st_mtime_ns, index_sha256=sha(Path(str(bam_path)+'.bai')),
        counts=dict(result), read_filters=dict(stats), templates_per_fragment_geometry=dict(family_sizes),
        phase_information_fraction_any_callable=n/denominator if denominator else None,
        direct_cocoverage_fraction_geometric_span_among_fetched=n/spanning if spanning else None,
        disposition='direct_molecular_cocoverage_observed_no_genotype_or_interaction_claim' if n else 'no_callable_direct_span_phase_unresolved',
        read_quality_orientation_diagnostics=[dict(site=i,allele=a,strand=s,read_arm=r,MAPQ=m,baseQ=b,reads=count)
            for (i,a,s,r,m,b),count in sorted(read_diagnostics.items())])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise SystemExit('Run only in authorized compute job')
    args.out.mkdir(parents=True, exist_ok=False, mode=0o700)
    shutil.copy2(__file__, args.out/'executed_data_read_phase.py')
    import pysam
    test = self_tests(pysam)
    manifest = source_manifest()
    table(args.out/'fixed_candidate_manifest.tsv', manifest)
    results = []
    with pysam.FastaFile(str(FASTA)) as fasta:
        for item in manifest:
            result = scan(item, fasta, pysam)
            write_json(args.out/(item['source']+'_'+item['donor']+'.json'), result)
            results.append(result)
    summary = dict(status='targeted_ATAC_phase_information_census_complete', job_id=os.environ['SLURM_JOB_ID'],
        checks=[test,'exact_historical_pair_donor_and_allelic_UID_join','reference_alleles_and_BAM_contig_lengths'],
        MAPQ_minimum=MAPQ, baseQ_minimum=BASEQ, BAMs=len(results), source_candidate_pairs=3,
        results=[{k:r[k] for k in ['identity','BAM_sample_identity_state','counts','phase_information_fraction_any_callable','direct_cocoverage_fraction_geometric_span_among_fetched','disposition']} for r in results],
        deduplication='Within each donor, combine primary proper-pair mates sharing RG/query_name; then collapse identical barcode/chromosome/outer-template coordinates across query names. Conflicting calls are excluded, not majority-voted. Same-coordinate barcode collisions can undercount separate molecules.',
        co_coverage='Both alleles must be called within the same read template before molecule deduplication. Geometric spanning alone, different templates, and deleted/skipped/low-quality bases cannot supply a phase observation.',
        interval_denominator='Every denominator is conditional on at least one alignment returned by indexed fetch over [first_variant0,last_variant0+1). Fragments with both sequenced mates outside this interval are not fetched even if their unsequenced insert spans both targets; geometric-span counts are not an all-spanning-fragment census.',
        barcode_scope='Cell Ranger corrected CB tag for GSE281367; source-prepended raw barcode prefix for GSE244832. Raw-barcode errors and source cell-QC membership are not newly resolved; this is donor-wide BAM information availability, not a lineage-specific analysis.',
        duplicate_limit='Custom GSE244832 BAMs did not establish PCR duplicate marking; geometry/barcode collapse is applied explicitly. Across-library barcode reuse remains an unverified collision risk.',
        mapping_bias='Quality, orientation, third-allele and allele-balance summaries are diagnostics, not allele-swapped remapping or a mapping-bias correction.',
        inferential_unit='Six prespecified donor-by-pair information summaries; fragments are not independent donors.',
        limitations='Selected historical candidate pairs and read-ascertained heterozygosity; no independent genotype truth, no disease response, no validation of ASE using the same reads, no interaction test, no GNMT reinterpretation.',
        genotype_download=False, raw_download=False, protected_outcomes_read=False, tests_of_interaction=False,
        source_sha256={str(p):sha(p) for p in [F1,Path(__file__),Path(str(FASTA)+'.fai'),
            ROOT/'data/GSE281367/metadata/donor_pairing.csv', ROOT/'data/GSE244832/metadata/donor_pairing.csv',
            ASE/'allelic_donor_counts.tsv.gz', ASE/'gse244832_allelic_donor_counts.tsv.gz']},
        environment=dict(python=platform.python_version(),pysam=pysam.__version__,node=platform.node()))
    write_json(args.out/'receipt.json',summary)
    size=sum(p.stat().st_size for p in args.out.iterdir() if p.is_file())
    if size>500000000: raise ValueError('Output cap exceeded; preserve diagnostics')
    print(json.dumps(dict(status=summary['status'],bytes=size,BAMs=len(results))),flush=True)


if __name__ == '__main__': main()
