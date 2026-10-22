# Figure guidance index

The authoritative project-wide figure rules remain in
[../FIGURE_GUIDELINES.md](../FIGURE_GUIDELINES.md). They cover panel density,
color semantics, typography, biological-unit display, PDF output, source-table
requirements, and archive conventions.

Additional authorities:

| Need | Location |
|---|---|
| Figure story and role | [../PAPER.md](../PAPER.md) |
| Figure readiness | [../STATUS.md](../STATUS.md) |
| Exact current values | [../RESULTS.md](../RESULTS.md) |
| Per-figure source inventory | `figures/README.md` and each figure directory’s `README.md` |
| Final reproducibility checks | [../ROADMAP.md](../ROADMAP.md) |

Every final plot must be generated from a read-only source table and reproduce
its values exactly. Control is gray (`#9E9E9E`), output is PDF with
`useDingbats=FALSE`, and cells/spots/runs must never visually imply a larger
biological sample size than the donor-level analysis.
