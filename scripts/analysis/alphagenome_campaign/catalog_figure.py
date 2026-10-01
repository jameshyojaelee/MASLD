#!/usr/bin/env python3
"""Candidate Figure 6 query specimen from the actual read-only evidence index."""
import argparse
import csv
import json
from pathlib import Path
import sqlite3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from catalog_query import query


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    with sqlite3.connect(a.db.resolve().as_uri()+'?mode=ro', uri=True) as db:
        rows = []
        for kind in ('variant', 'gene', 'region', 'event', 'program', 'context'):
            identity = db.execute('SELECT identifier FROM entity WHERE kind=? ORDER BY identifier LIMIT 1', (kind,)).fetchone()[0]
            result = query(a.db, kind, identity, limit=10)
            (a.out/(kind+'.json')).write_text(json.dumps(result, indent=2)+'\n')
            n = result['total_evidence']
            rows.append(dict(entry=kind, identity=identity, evidence_records=n,
                             shown_records=len(result['evidence']), relationships=len(result['relationships']),
                             status=result['status']))
    with (a.out/'Figure6_query_values.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':6, 'pdf.fonttype':42,
                         'axes.linewidth':0.5, 'savefig.transparent':False})
    fig, ax = plt.subplots(figsize=(7.1, 2.3))
    ax.set_axis_off()
    cell_rows = [[r['entry'].capitalize(), r['identity'], str(r['evidence_records']), str(r['relationships'])] for r in rows]
    table = ax.table(cellText=cell_rows, colLabels=['Entry', 'Example identity', 'Evidence records', 'Relationships'],
                     cellLoc='left', colLoc='left', colWidths=[.10,.64,.14,.12], bbox=[0,.24,1,.64])
    table.auto_set_font_size(False); table.set_fontsize(6)
    for (i,j), cell in table.get_celld().items():
        cell.set_linewidth(0.3); cell.set_edgecolor('#D9D9D9')
        cell.set_facecolor('#F0F0F0' if i == 0 else 'white')
        cell.PAD=.035
    ax.text(0,.98,'Candidate Figure 6 | Six molecular query entry points', transform=ax.transAxes, va='top', fontsize=6)
    ax.text(0,.17,'Examples selected by smallest stable identity, without effect-size selection. Counts include recorded and unsupported evidence.',
            transform=ax.transAxes, va='top', fontsize=6)
    ax.text(0,.10,'Program answers retain fixed membership links; no program-level molecular effects are inferred. Full context totals exceed the 10-record example display.',
            transform=ax.transAxes, va='top', fontsize=6)
    ax.text(0,.03,'Local aggregate attachment specimen; measured and predicted effects retain separate units and source conditions. Candidate, not adopted.',
            transform=ax.transAxes, va='top', fontsize=6)
    fig.subplots_adjust(left=.03,right=.99,top=.99,bottom=.03)
    fig.savefig(a.out/'Figure6_query_specimen.pdf')
    plt.close(fig)


if __name__ == '__main__':
    main()
