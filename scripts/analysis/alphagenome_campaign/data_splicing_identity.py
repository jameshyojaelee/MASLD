#!/usr/bin/env python3
"""Source-native splice identity/readiness, without new outcome selection or inference.

The available significant-pair deposit is an outcome-selected SOURCE population.
This script admits all its rows and applies identity-only readiness rules. It
does not reconstruct the full tested event population or full LeafCutter clusters.
"""
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
import shutil

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data"
SOURCE = ROOT / "data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL"
OLD = ROOT / "GWAS/finemapping/results/alphagenome_atlas/p3a-splice-direction-20260915T141415Z/tables"
PAIRS = SOURCE / "Liver.v8.sqtl_signifpairs.txt.gz"
SGENES = SOURCE / "Liver.v8.sgenes.txt.gz"


def rows(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def source_event(identity):
    chrom, start, end, cluster, gene = identity.split(':')
    start, end = int(start), int(end)
    if start < 0 or end <= start or not gene.startswith('ENSG'):
        raise ValueError('Invalid source event identity: ' + identity)
    return chrom, start, end, cluster, gene


def resolve_strand(identity, by_event, by_gene):
    chrom, _, _, _, gene = source_event(identity)
    values = by_event.get(identity, set())
    state = 'exact_source_phenotype'
    if not values:
        values = by_gene.get(gene, set())
        state = 'exact_versioned_source_gene'
    if len(values) != 1:
        return '.', 'missing_annotation' if not values else 'conflicting_source_annotations'
    gchrom, strand = next(iter(values))
    if gchrom != chrom or strand not in ('+', '-'):
        return '.', 'chromosome_or_strand_conflict'
    return strand, state


def duplicate_details(records):
    fields = set().union(*(r.keys() for r in records))
    differs = sorted(k for k in fields if len({r.get(k, '') for r in records}) > 1)
    return len(records), 'exact_repeated_record' if not differs else 'conflicting_records', ';'.join(differs)


def self_tests():
    event = 'chr1:10:20:clu_1:ENSG00000000001.2'
    assert resolve_strand(event, {event: {('chr1', '-')}}, {}) == ('-', 'exact_source_phenotype')
    assert resolve_strand(event, {}, {'ENSG00000000001.2': {('chr1', '+')}}) == ('+', 'exact_versioned_source_gene')
    assert resolve_strand(event, {}, {'ENSG00000000001.3': {('chr1', '+')}})[0] == '.'
    assert resolve_strand(event, {}, {'ENSG00000000001.2': {('chr1', '+'), ('chr1', '-')}})[0] == '.'
    assert resolve_strand(event, {}, {'ENSG00000000001.2': {('chr2', '+')}})[0] == '.'
    a = {'state': 'scored', 'slope': '0.5', 'effect': '0.1'}
    assert duplicate_details([a, dict(a)]) == (2, 'exact_repeated_record', '')
    assert duplicate_details([a, dict(a, slope='0.6')]) == (2, 'conflicting_records', 'slope')
    assert duplicate_details([a, dict(a, state='missing')]) == (2, 'conflicting_records', 'state')
    return 'exact_gene_version_no_nearest_gene_no_conflicting_strand_and_duplicate_record_cases_passed'


def table(path, records, columns):
    with path.open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        w.writeheader()
        w.writerows(records)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Run identity checks on a compute node through the authorized launcher')
    args.out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.out / 'executed_data_splicing_identity.py')
    checks = [self_tests()]
    by_event, by_gene = defaultdict(set), defaultdict(set)
    for r in rows(SGENES):
        parsed = source_event(r['phenotype_id'])
        if parsed[4] != r['gene_id']:
            raise ValueError('Source sgenes phenotype and gene_id disagree')
        value = (r['gene_chr'], r['strand'])
        by_event[r['phenotype_id']].add(value)
        by_gene[r['gene_id']].add(value)
    legacy = defaultdict(list)
    legacy_indices = defaultdict(list)
    incomplete_legacy = []
    legacy_counts = Counter()
    for index, r in enumerate(rows(OLD / 'splice_direction_pairs.tsv')):
        legacy_counts[r['state']] += 1
        if r.get('phenotype_id'):
            identity = (r['variant_uid'], r['phenotype_id'])
            legacy[identity].append(r)
            legacy_indices[identity].append(index)
        else:
            incomplete_legacy.append(dict(historical_row=index, variant_uid=r['variant_uid'], state=r['state'],
                missing='source_phenotype_id_absent_in_historical_failed_row;do_not_assign_an_event_by_proximity'))
    table(args.out / 'historical_rows_without_event_identity.tsv', incomplete_legacy,
          ['historical_row', 'variant_uid', 'state', 'missing'])
    repeated_legacy = []
    for (uid, phenotype), records in legacy.items():
        if len(records) > 1:
            multiplicity, state, differs = duplicate_details(records)
            repeated_legacy.append(dict(variant_uid=uid, source_phenotype_id=phenotype,
                historical_row_indices=';'.join(map(str, legacy_indices[(uid, phenotype)])),
                multiplicity=multiplicity, duplicate_state=state, differing_columns=differs))
    table(args.out / 'historical_duplicate_identities.tsv', repeated_legacy,
        ['variant_uid', 'source_phenotype_id', 'historical_row_indices', 'multiplicity', 'duplicate_state', 'differing_columns'])
    source_counts, source_fingerprints = Counter(), defaultdict(set)
    for r in rows(PAIRS):
        key = (r['variant_id'], r['phenotype_id'])
        source_counts[key] += 1
        # Full-record equality is only a duplicate-evidence check; no effect
        # value selects a row or changes its scientific readiness.
        source_fingerprints[key].add(hashlib.sha256(json.dumps(r, sort_keys=True).encode()).hexdigest())
    repeated_source = [dict(source_variant_id=k[0], source_phenotype_id=k[1], multiplicity=n,
        duplicate_state='exact_repeated_record' if len(source_fingerprints[k]) == 1 else 'conflicting_records')
        for k, n in source_counts.items() if n > 1]
    table(args.out / 'source_duplicate_identities.tsv', repeated_source,
        ['source_variant_id', 'source_phenotype_id', 'multiplicity', 'duplicate_state'])
    fields = ['source_row_index', 'source_variant_id', 'source_phenotype_id', 'source_identity_multiplicity',
        'variant_chrom', 'variant_position1', 'ref', 'alt',
        'event_chrom', 'intron_start_source', 'intron_end_source', 'source_cluster_id', 'gene_id_versioned',
        'source_gene_strand', 'strand_identity_source', 'allele_type', 'legacy_state', 'legacy_identity_multiplicity',
        'source_cluster_complete', 'raw_REF_ALT_junctions_available', 'new_effect_evaluation_ready', 'missing_codes']
    counts, strand_states, allele_types = Counter(), Counter(), Counter()
    cluster_members, cluster_strands = defaultdict(set), defaultdict(set)
    seen_pairs, seen_legacy = set(), set()
    events = {}
    manifest = args.out / 'splicing_event_manifest_v1.tsv'
    with manifest.open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for source_index, r in enumerate(rows(PAIRS)):
            variant_id, phenotype = r['variant_id'], r['phenotype_id']
            vchrom, vpos, ref, alt, assembly = variant_id.split('_')
            if assembly != 'b38' or int(vpos) < 1:
                raise ValueError('Source variant build/coordinate differs')
            chrom, start, end, cluster, gene = source_event(phenotype)
            strand, state = resolve_strand(phenotype, by_event, by_gene)
            uid = f'{vchrom}:{vpos}:{ref}:{alt}'
            key = (variant_id, phenotype)
            seen_pairs.add(key)
            allele_type = 'SNV' if len(ref) == len(alt) == 1 and ref != alt and set(ref + alt) <= set('ACGT') else 'indel_or_other'
            codes = ['C', 'J', 'R']
            if strand == '.': codes.append('S')
            if allele_type != 'SNV': codes.append('A')
            if vchrom != chrom: codes.append('X')
            if len(source_fingerprints[key]) > 1: codes.append('D')
            oldkey = (uid, phenotype)
            oldrecords = legacy.get(oldkey, [])
            oldstates = sorted({x['state'] for x in oldrecords})
            oldstate = ';'.join(oldstates) if oldstates else 'not_in_historical_targeted_rows'
            if oldkey in legacy: seen_legacy.add(oldkey)
            w.writerow(dict(source_row_index=source_index, source_variant_id=variant_id, source_phenotype_id=phenotype,
                source_identity_multiplicity=source_counts[key],
                variant_chrom=vchrom, variant_position1=vpos, ref=ref, alt=alt, event_chrom=chrom,
                intron_start_source=start, intron_end_source=end, source_cluster_id=cluster,
                gene_id_versioned=gene, source_gene_strand=strand, strand_identity_source=state,
                allele_type=allele_type, legacy_state=oldstate, legacy_identity_multiplicity=len(oldrecords), source_cluster_complete=False,
                raw_REF_ALT_junctions_available=False, new_effect_evaluation_ready=False,
                missing_codes=';'.join(codes)))
            counts['source_variant_event_pairs'] += 1
            strand_states[state] += 1
            allele_types[allele_type] += 1
            ckey = (chrom, cluster)
            cluster_members[ckey].add((start, end))
            if strand != '.': cluster_strands[ckey].add(strand)
            events[phenotype] = dict(source_phenotype_id=phenotype, event_chrom=chrom,
                intron_start_source=start, intron_end_source=end, source_cluster_id=cluster,
                gene_id_versioned=gene, source_gene_strand=strand, strand_identity_source=state)
    if legacy.keys() - seen_legacy:
        raise ValueError('Historical exact event identities absent from source: ' + str(len(legacy.keys() - seen_legacy)))
    checks.append('all_historical_rows_with_exact_event_identity_join_source')
    table(args.out / 'unique_source_events.tsv', events.values(), list(next(iter(events.values()))))
    cluster_rows = [dict(event_chrom=c, source_cluster_id=k, observed_significant_source_introns=len(v),
        gene_annotation_strands=';'.join(sorted(cluster_strands[(c, k)])), complete_source_denominator=False)
        for (c, k), v in sorted(cluster_members.items())]
    table(args.out / 'observed_cluster_membership.tsv', cluster_rows, list(cluster_rows[0]))
    missing = dict(
        C='Complete GTEx v8 liver LeafCutter cluster membership including nonsignificant/filtered-out-of-QTL-test members, original clustering convention and source intron coordinate convention; observed significant members are not a complete denominator.',
        J='New raw per-arm splice-junction coordinates AND strand, track IDs/order, per-track usage, window/model identity, and original proposed-junction lists; preserve REF union ALT and explicit absent/present flags. Archived aggregate counts cannot reconstruct gained/lost junctions.',
        R='Reference-allele validation and exact coordinate-conversion probe for the new inference route. Historical model junction mapping was source(start,end) -> (start,end-1), measured in a specific API version; do not silently assume every route shares it.',
        S='Unique exact source phenotype or versioned-gene strand annotation absent/conflicting. Never substitute nearest-gene strand.',
        A='Initial bounded scalar splice comparison supports SNVs only; indels/other alleles need normalization and a separate sequence-to-junction coordinate mapping.',
        X='Variant and event chromosomes differ; verify source identity before admission.',
        D='Repeated exact source variant/event identity has conflicting source records; retain every row, do not select one using outcomes.')
    receipt = dict(status='partial_source_identity_manifest_exact_effect_comparison_not_ready', job_id=os.environ['SLURM_JOB_ID'],
        node=platform.node(), checks=checks, counts=dict(counts), unique_source_events=len(events),
        observed_source_clusters=len(cluster_members), strand_join_rows=dict(strand_states), allele_types=dict(allele_types),
        clusters_with_conflicting_gene_annotation_strands=sum(len(s) > 1 for s in cluster_strands.values()),
        historical_table_rows=sum(legacy_counts.values()), historical_states=dict(legacy_counts),
        historical_rows_with_event_identity=sum(len(v) for v in legacy.values()),
        historical_distinct_variant_event_identities=len(legacy), historical_failed_rows_without_event_identity=len(incomplete_legacy),
        historical_duplicate_identities=len(repeated_legacy),
        historical_extra_duplicate_rows=sum(len(v)-1 for v in legacy.values()),
        historical_duplicate_record_states=dict(Counter(r['duplicate_state'] for r in repeated_legacy)),
        historical_distinct_scored_identities=sum(any(r['state'] == 'scored' for r in v) for v in legacy.values()),
        source_distinct_variant_event_identities=len(seen_pairs), source_duplicate_identities=len(repeated_source),
        source_extra_duplicate_rows=sum(n-1 for n in source_counts.values()),
        source_duplicate_record_states=dict(Counter(r['duplicate_state'] for r in repeated_source)),
        historical_distinct_exact_identities_joined=len(seen_legacy), new_effect_evaluation_ready_rows=0,
        eligibility='No p-value, slope, effect size or prediction value used for new row selection; every available source significant pair retained. The upstream significant-pair population is outcome-selected and is not the all-tested universe.',
        strand_scope='Exact source-gene strand annotation, with exact phenotype or exact versioned-gene linkage; not an independently stranded junction assay. Full source cluster direction still requires source clustering definition.',
        old_comparator='Historical direct start/end neighbor denominator on REF-intersection-ALT predicted junctions; coordinates omit strand. Preserved as loose comparator, not exact source cluster; no rerun.',
        gained_lost_future_definition='Join raw REF and ALT proposed junctions on chromosome,start,end,strand,track; preserve full union and availability flags. Define gained/lost as arm membership, not evidence of biological zero usage when unreturned; fixed source-cluster contrast and novel-junction sensitivity remain distinct.',
        native_effect_units='GTEx normalized intron-excision phenotype per ALT dosage; not deltaPSI; no effect values copied to readiness manifest',
        missing_field_codes=missing, new_inference_run=False, protected_outcomes_read=False,
        input_sha256={str(p):sha(p) for p in [PAIRS, SGENES, OLD/'splice_direction_pairs.tsv', OLD/'splice_direction.json', Path(__file__)]})
    (args.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    output_bytes = sum(p.stat().st_size for p in args.out.iterdir() if p.is_file())
    if output_bytes > 100_000_000:
        raise ValueError('Output exceeds authorized .1GB cap; preserve files and report without deleting')
    print(json.dumps(dict(receipt=receipt, output_bytes=output_bytes), indent=2), flush=True)


if __name__ == '__main__':
    main()
