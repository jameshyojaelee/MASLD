#!/usr/bin/env python3
"""Prepare the approved Figure S5P–AA additions from completed, fixed analyses."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import render_final_candidate as base
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'Render on a compute node.'
    c = args.candidate.resolve()
    context = c / 'disease_context'
    assert json.loads((c / 'final_review.json').read_text())['status'] == 'completed'
    out = args.output.resolve()
    out.mkdir(exist_ok=False)
    panels = out / 'panels'
    sys.argv = [__file__, '--candidate', str(c), '--context', str(context), '--output', str(panels)]
    base.main()
    data = out / 'data'
    data.mkdir()
    saved = []

    def export(source, name, interpretation):
        source = source.resolve()
        dest = data / name
        assert not dest.exists()
        shutil.copy2(source, dest)
        assert digest(source) == digest(dest)
        saved.append(dict(file=name, source=str(source), sha256=digest(dest), interpretation=interpretation))

    # KEY MESSAGE: batch adjustment changes discovery counts on exactly matched donors/peaks.
    counts = pd.read_csv(c / 'batch_comparison.tsv', sep='\t')
    expected = {'hepatocyte': (1, 3), 'stellate': (820, 291)}
    for row in counts.itertuples():
        assert (row.n_condition_BH_q05, row.n_batch_BH_q05) == expected[row.lineage]
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    effects = pd.read_csv(context / 'accessibility_effects.tsv.gz', sep='\t')
    for lineage, color, marker in [('hepatocyte', base.BLUE, 'o'), ('stellate', base.MAGENTA, '^')]:
        g = effects[effects.cohort.eq('GSE244832') & effects.lineage.eq(lineage) & effects.batch_logFC.notna()]
        ax.scatter(g.logFC, g.batch_logFC, s=1, alpha=.2, c=color, marker=marker,
                   rasterized=True, label=lineage)
    lo, hi = ax.get_xlim()
    ax.plot([lo, hi], [lo, hi], c=base.GRAY, lw=.6)
    ax.axhline(0, c='black', lw=.4)
    ax.axvline(0, c='black', lw=.4)
    ax.text(.02, .97, 'BH calls: condition-only → batch-adjusted\nHepatocyte: 1 → 3; stellate: 820 → 291',
            transform=ax.transAxes, va='top')
    ax.set_xlabel('Condition-only log2 fragment effect')
    ax.set_ylabel('Recorded-batch-adjusted effect')
    ax.legend(frameon=False, markerscale=4, loc='lower right')
    fig.tight_layout()
    fig.savefig(panels / '04_batch_sensitivity.pdf', bbox_inches='tight')
    plt.close(fig)

    # KEY MESSAGE: gene names in examples identify descriptive links, not regulated targets.
    examples = pd.read_csv(context / 'deterministic_examples.tsv', sep='\t')
    labels = {'peakwise_candidate': 'Peakwise candidate', 'source_dependent': 'Source-dependent',
              'indeterminate': 'Indeterminate', 'batch_sensitive': 'Batch-sensitive'}
    fig, ax = plt.subplots(figsize=(4.6, 2.6))
    ax.errorbar(examples.logFC, np.arange(len(examples)),
                xerr=[examples.logFC - examples.ci_lower, examples.ci_upper - examples.logFC],
                fmt='o', ms=3, c=base.MAGENTA, lw=.7)
    ax.axvline(0, c=base.GRAY, lw=.6)
    ax.set_yticks(range(len(examples)), [f'{labels[r.display_category]}; link: {r.gene_symbol}\n{r.peak_coordinate}'
                                      for r in examples.itertuples()])
    ax.set_xlabel('MASH vs normal log2 fragment effect (95% moderated-t CI)')
    fig.tight_layout()
    fig.savefig(panels / '08_example_uncertainty.pdf', bbox_inches='tight')
    plt.close(fig)

    # KEY MESSAGE: small internal predictive gains persist through the bounded stress tests.
    comparison = pd.read_csv(c / 'prediction_summary/comparisons.tsv', sep='\t')
    names = list(dict.fromkeys(comparison.scenario))
    display = {'baseline': 'Original wider grid', 'cap_500': '500 nuclei per donor'}
    display.update({f'aggregation_{s}': f'Aggregation seed {s}' for s in (43, 44)})
    display.update({f'pairing_shuffle_{s}': f'RNA-pairing shuffle {s}' for s in (43, 44, 45)})
    display.update({f'omit_training_{s}': f'Training omission: donor {s}' for s in (342, 786, 470)})
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for method, color, offset, marker in [('promoter_only', base.GRAY, -.1, 'o'), ('distance', base.MAGENTA, .1, 's')]:
        g = comparison[comparison.baseline.eq(method)].set_index('scenario').loc[names]
        ax.errorbar(g.mean_gain, np.arange(len(names)) + offset,
                    xerr=[g.mean_gain - g.ci_low, g.ci_high - g.mean_gain], fmt=marker, ms=3,
                    c=color, lw=.6, label='Promoter-only' if method == 'promoter_only' else 'Distance')
    ax.set_yticks(range(len(names)), [display[s] for s in names])
    ax.axvline(0, c='black', lw=.5)
    ax.set_xlabel('Baseline minus Cicero standardized MSE\n(95% conditional donor interval; 39 held-out donors)')
    ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(.5, 1.10), ncol=2)
    fig.tight_layout()
    fig.savefig(panels / '09_prediction_robustness.pdf', bbox_inches='tight')
    plt.close(fig)

    # KEY MESSAGE: correct pairing has no demonstrated Cicero advantage beyond proximity.
    changes = pd.read_csv(c / 'prediction_summary/paired_scenario_change_summary.tsv', sep='\t')
    g = changes[changes.scenario.str.startswith('pairing_shuffle') & changes.baseline.eq('distance')]
    assert len(g) == 3 and (g.ci_low < 0).all() and (g.ci_high > 0).all()
    fig, ax = plt.subplots(figsize=(3.9, 2.3))
    ax.errorbar(g.original_gain_minus_scenario_gain, range(len(g)),
                xerr=[g.original_gain_minus_scenario_gain - g.ci_low,
                      g.ci_high - g.original_gain_minus_scenario_gain], fmt='o', ms=3, c=base.MAGENTA, lw=.7)
    ax.set_yticks(range(len(g)), [display[s] for s in g.scenario])
    ax.axvline(0, c=base.GRAY, lw=.6)
    ax.set_xlabel('Paired minus shuffled Cicero advantage over distance\n(standardized MSE; 95% conditional donor interval)')
    fig.tight_layout()
    fig.savefig(panels / '10_pairing_controls.pdf', bbox_inches='tight')
    plt.close(fig)

    # KEY MESSAGE: two hepatocyte RNA candidates pass both matched models; batch can introduce calls.
    rna = pd.read_csv(c / 'rna_models/rna_effects.tsv.gz', sep='\t')
    genes = ['CHST9', 'TMEM150B', 'TRIM9']
    selected = rna[rna.gene_symbol.isin(genes)].copy()
    assert len(selected) == 6
    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    for adjustment, color, offset, marker in [('condition_only', base.GRAY, -.13, 'o'), ('batch_jb', base.MAGENTA, .13, 's')]:
        g = selected[selected.adjustment.eq(adjustment)].set_index('gene_symbol').loc[genes]
        y = np.arange(3) + offset
        for i, (gene, row) in enumerate(g.iterrows()):
            ax.errorbar(row.logFC, y[i], xerr=[[row.logFC - row.ci_lower], [row.ci_upper - row.logFC]],
                        fmt=marker, ms=3.5, color=color, markerfacecolor=color if row.qvalue < .05 else 'white', lw=.7)
            ax.text(5.7, y[i], f'q={row.qvalue:.4f}', color=color, va='center')
    ax.axvline(0, c=base.GRAY, lw=.6)
    ax.set_yticks(range(3), genes)
    ax.set_xlim(-.15, 8.0)
    ax.invert_yaxis()
    ax.set_xlabel('MASH vs normal log2 normalized RNA UMI effect\n(pointwise 95% moderated-t CI; 5 normal / 9 MASH donors)')
    handles = [Line2D([], [], color=base.GRAY, marker='o', ls='', label='Condition-only'),
               Line2D([], [], color=base.MAGENTA, marker='s', ls='', label='Recorded batch')]
    ax.legend(handles=handles, frameon=False, loc='upper center', bbox_to_anchor=(.5, 1.22), ncol=2)
    ax.text(.02, 1.03, 'Filled: BH q < 0.05 across 21,054 tested genes', transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig(panels / '11_RNA_disease_effects.pdf', bbox_inches='tight')
    plt.close(fig)
    selected.to_csv(data / 'displayed_RNA_candidates.tsv', sep='\t', index=False)

    peak = pd.read_csv(context / 'cross_cohort_peak_evidence.tsv.gz', sep='\t')
    quartet = peak[peak.peak_evidence_state.eq('two_cohort_peakwise_candidate')]
    assert len(quartet) == 4 and quartet.lineage.eq('stellate').all()
    quartet.to_csv(data / 'four_stellate_peakwise_candidates.tsv', sep='\t', index=False)
    export(context / 'accessibility_effects.tsv.gz', 'accessibility_effects.tsv.gz', 'Peak effects; peak-level BH within cohort, lineage and model.')
    for name in ['cohort_agreement.tsv', 'cross_cohort_peak_evidence.tsv.gz', 'cross_cohort_link_agreement.tsv',
                 'deterministic_examples.tsv', 'fixed_program_link_coverage.tsv',
                 'fixed_program_membership_annotation.tsv.gz', 'gene_link_coverage.tsv.gz', 'peak_to_gene_context.tsv.gz']:
        export(context / name, name, 'Descriptive accessibility context; link scores are not probabilities; no gene-level regulatory call.')
    for name in ['batch_comparison.tsv', 'rna_recount_checks.tsv', 'rna_molecule_accounting.tsv']:
        export(c / name, name, 'Matched sensitivity or exact reconstruction check.')
    for part, names_to_copy in {
        'disease_synthesis': ['roster.tsv', 'influence.tsv', 'models.tsv', 'overlaps.tsv'],
        'rna_pseudobulk': ['donors.tsv', 'source_diagnosis_decisions.tsv', 'library_decisions.tsv'],
        'rna_models': ['rna_effects.tsv.gz', 'model_decisions.tsv'],
        'donor_associations': ['decisions.tsv'],
        'prediction_summary': ['comparisons.tsv', 'donor_scores.tsv', 'coverage.tsv', 'scenario_key_intersections.tsv',
                               'paired_scenario_change_summary.tsv', 'tuning_boundaries.tsv'],
    }.items():
        for name in names_to_copy:
            export(c / part / name, f'{part}_{name}', 'Donor-level analysis; original and sensitivity hypothesis families remain separate.')
    for source in sorted((c / 'prediction_summary').glob('*excluded_genes.tsv.gz')):
        export(source, 'prediction_' + source.name, 'Explicit excluded gene/fold IDs; no mechanism inferred from absent scores.')
    for source in sorted((c / 'motif_context_index_corrected').glob('*/annotation/motif_family_gene.tsv.gz')):
        export(source, source.parents[1].name + '_motif_family_gene.tsv.gz', 'Motif-family compatibility only, not identified TF activity or binding.')
    earlier = c.parents[1] / '20261024'
    for part, names_to_copy in {
        'ridge_grid_sensitivity_r1/summary': ['comparisons.tsv', 'coverage.tsv', 'grid_effects.tsv', 'original_grid_reference.tsv', 'ridge_alpha_distribution.tsv'],
        'wide_proximity_comparison_r1': ['comparisons.tsv', 'fold_coverage.tsv', 'common_genes.tsv'],
    }.items():
        for name in names_to_copy:
            export(earlier / part / name, part.replace('/', '_') + '_' + name,
                   'Earlier completed comparison; retain its original multiplicity correction separately from new descriptive stress tests.')
    pd.DataFrame(saved).to_csv(data / 'source_tables.tsv', sep='\t', index=False)

    index = pd.read_csv(panels / 'panels.tsv', sep='\t')
    letters = list('pqrstuvwxyz') + ['aa']
    assert len(index) == len(letters) == 12
    index['callout'] = ['S5' + s.upper() for s in letters]
    index['manuscript_file'] = ['figS5' + s + '.pdf' for s in letters]
    index['sha256'] = [digest(panels / n) for n in index.file]
    index.to_csv(out / 'panel_index.tsv', sep='\t', index=False)
    checks = dict(status='completed', n_panels=12, four_stellate_candidates=4, displayed_RNA_rows=6,
                  fixed_primary_method='TMM_voom_quality_weights', program_definitions_changed=False,
                  catalog_classes_changed=False, new_model_fits=False, slurm_job_id=os.environ['SLURM_JOB_ID'])
    (out / 'completion.json').write_text(json.dumps(checks, indent=2))
    print(json.dumps(checks))


if __name__ == '__main__':
    main()
