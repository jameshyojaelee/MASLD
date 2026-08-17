#!/usr/bin/env python
"""T1 helper: exact cells-per-run-per-cell-type from the Hotspot cell-score tables.

WHY THIS EXISTS. No donor-level program-score table in this project carries the
number of cells the score was averaged over, and a donor program score computed
from 12 cells is not the same measurement as one computed from 3,000. The only
substrate that gives the exact denominator actually used by the scoring step is
the per-cell score table itself, keyed on "{SRR}_{barcode}". The
phase05 lineage count table is a different annotation pass over the atlas and
disagrees with the scored cell set (it reports 0 T cells for donors that carry a
T-cell program score), so it is used only for the abundance denominator, never
for the precision gate.

Emits one row per (cell_type, run, module-invariant cell count). Run-to-donor
collapse happens in R with the shared lib_donor_collapse.R pairing tables.
"""

import sys
import pyarrow.parquet as pq
import pyarrow.compute as pc

CELL_TYPES = ["cholangiocytes", "fibroblasts", "hepatocytes", "macrophages", "tcells"]


def main(hotspot_root: str, out_path: str) -> None:
    rows = []
    for cell_type in CELL_TYPES:
        path = f"{hotspot_root}/{cell_type}/cell_scores.parquet"
        table = pq.read_table(path, columns=["cell_id", "module"])
        modules = pc.unique(table.column("module"))
        first = modules[0]
        one = table.filter(pc.equal(table.column("module"), first))
        cell_ids = one.column("cell_id").to_pylist()
        # A cell id is "{run}_{barcode}"; the barcode is the final underscore field.
        counts = {}
        for cid in cell_ids:
            run = cid.rsplit("_", 1)[0]
            counts[run] = counts.get(run, 0) + 1
        # The cell set must not depend on which module was chosen.
        n_per_module = len(one)
        for other in modules.to_pylist()[1:]:
            if len(table.filter(pc.equal(table.column("module"), other))) != n_per_module:
                raise SystemExit(
                    f"{cell_type}: module {other} scores a different cell set than {first}"
                )
        for run, n in sorted(counts.items()):
            rows.append((cell_type, run, n))
        print(f"{cell_type}\truns={len(counts)}\tcells={n_per_module}", file=sys.stderr)

    with open(out_path, "w") as handle:
        handle.write("cell_type\trun\tn_cells\n")
        for cell_type, run, n in rows:
            handle.write(f"{cell_type}\t{run}\t{n}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
