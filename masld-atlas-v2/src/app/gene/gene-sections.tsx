"use client";

import { type ReactNode } from "react";
import {
  ForestPlot,
  Heatmap,
  DotPlot,
  DivergingBar,
  type ForestRow,
  type HeatmapCell,
  type DotPlotPoint,
  type DivergingBarRow,
} from "@/components/charts";
import { DataTable, type DataTableColumn } from "@/components/data-table";
import { StatTile } from "@/components/stat-tile";
import { AncestryChip, StageChip, ProgramChip } from "@/components/chips";
import { EmptyState, SkeletonBlock } from "@/components/states";
import { sqlString } from "@/lib/duck";
import { devStageColor } from "@/lib/palette";
import {
  useParquetRows,
  type QueryResult,
  num,
  clamp01,
  fmtLogFC,
  fmtP,
  fmtPP4,
} from "./use-gene-data";

// ---------------------------------------------------------------------------
// Shared shells
// ---------------------------------------------------------------------------

/** Card wrapper with a title + optional description. */
function Panel({
  title,
  description,
  children,
}: {
  title: ReactNode;
  description?: ReactNode;
  children: ReactNode;
}) {
  return (
    <section className="rounded-lg border border-border bg-card p-4 shadow-sm">
      <h3 className="text-sm font-semibold uppercase tracking-wider text-muted-foreground">
        {title}
      </h3>
      {description && (
        <p className="mt-0.5 mb-3 text-xs text-muted-foreground">{description}</p>
      )}
      {!description && <div className="mb-3" />}
      {children}
    </section>
  );
}

/** Render loading / error / empty / data uniformly for a query result. */
function QueryBlock<T>({
  q,
  empty,
  children,
}: {
  q: QueryResult<T>;
  empty: ReactNode;
  children: (rows: T[]) => ReactNode;
}) {
  if (q.error)
    return (
      <p className="text-sm text-destructive">Failed to load data: {q.error}</p>
    );
  if (q.loading) return <SkeletonBlock className="h-56" />;
  if (!q.data || q.data.length === 0)
    return <EmptyState title={empty} />;
  return <>{children(q.data)}</>;
}

const where = (col: string, symbol: string) =>
  `WHERE upper(${col}) = upper(${sqlString(symbol)})`;

/** Natural stage ordering (F0<F1<…, else by embedded number, else lexicographic). */
function stageSort(a: string, b: string): number {
  const na = a.match(/(\d+)/);
  const nb = b.match(/(\d+)/);
  if (na && nb) return Number(na[1]) - Number(nb[1]);
  return a.localeCompare(b);
}

function groupBy<T>(rows: T[], key: (r: T) => string): Map<string, T[]> {
  const m = new Map<string, T[]>();
  for (const r of rows) {
    const k = key(r);
    const arr = m.get(k) ?? [];
    arr.push(r);
    m.set(k, arr);
  }
  return m;
}

/** Entry with the most rows. */
function mostPopulated<T>(m: Map<string, T[]>): [string, T[]] {
  return [...m.entries()].sort((a, b) => b[1].length - a[1].length)[0];
}

// Retired / non-canonical contrast labels that must never be the default view.
const LEGACY_CONTRAST = /legacy|dream/i;

// ---------------------------------------------------------------------------
// Expression: per-cohort forest + stage trajectory + (conditional) sex
// ---------------------------------------------------------------------------

interface CohortRow {
  dataset: string;
  contrast: string;
  logfc: number;
  padj: number;
}
interface TrajRow {
  axis: string;
  reference: string;
  stage: string;
  logfc: number;
  padj: number;
  n_cohorts: number;
}
interface SexRow {
  sex_class: string;
  logfc_female: number;
  logfc_male: number;
  interaction_padj: number;
}

export function ExpressionSection({ symbol }: { symbol: string }) {
  const cohorts = useParquetRows<CohortRow>(
    "gene_per_cohort_de.parquet",
    `SELECT dataset, contrast, logfc, padj FROM gene_per_cohort_de ${where(
      "symbol",
      symbol
    )}`
  );
  const traj = useParquetRows<TrajRow>(
    "gene_trajectories.parquet",
    `SELECT axis, reference, stage, logfc, padj, n_cohorts FROM gene_trajectories ${where(
      "symbol",
      symbol
    )}`
  );
  const sex = useParquetRows<SexRow>(
    "sex_gene_classification.parquet",
    `SELECT sex_class, logfc_female, logfc_male, interaction_padj FROM sex_gene_classification ${where(
      "symbol",
      symbol
    )}`
  );

  return (
    <div className="grid gap-4 lg:grid-cols-2">
      <Panel
        title="Per-cohort differential expression"
        description="Disease-vs-control log2 fold change in each contributing cohort."
      >
        <QueryBlock q={cohorts} empty="No per-cohort estimates for this gene.">
          {(rows) => {
            // Exclude retired contrasts; prefer the canonical disease-vs-control,
            // else the most-populated remaining contrast.
            const usable = rows.filter((r) => !LEGACY_CONTRAST.test(r.contrast));
            const pool = usable.length ? usable : rows;
            const byContrast = groupBy(pool, (r) => r.contrast);
            const [contrast, picked] = byContrast.has("disease_vs_control")
              ? ["disease_vs_control", byContrast.get("disease_vs_control")!]
              : mostPopulated(byContrast);
            const data: ForestRow[] = picked
              .slice()
              .sort((a, b) => b.logfc - a.logfc)
              .map((r) => ({
                label: r.dataset,
                estimate: num(r.logfc) ?? 0,
                sig: (num(r.padj) ?? 1) < 0.05,
              }));
            return (
              <ForestPlot
                data={data}
                xLabel="log2 fold change"
                valueFormat={(v) => fmtLogFC(v)}
                caption={`Contrast: ${contrast}`}
                height={Math.max(160, 40 + data.length * 22)}
                ariaLabel={`Per-cohort fold change for ${symbol}`}
              />
            );
          }}
        </QueryBlock>
      </Panel>

      <Panel
        title="Stage trajectory"
        description="Fold change across disease stages relative to the reference."
      >
        <QueryBlock q={traj} empty="No stage-trajectory estimates for this gene.">
          {(rows) => {
            // Pick the axis with the most rows, then a single reference within
            // it (avoids duplicate stage bars from multiple baselines), then one
            // row per stage in natural order.
            const [axis, axisRows] = mostPopulated(groupBy(rows, (r) => r.axis));
            const [ref, refRows] = mostPopulated(
              groupBy(axisRows, (r) => r.reference)
            );
            const seen = new Set<string>();
            const data: DivergingBarRow[] = refRows
              .slice()
              .sort((a, b) => stageSort(a.stage, b.stage))
              .filter((r) => {
                if (seen.has(r.stage)) return false;
                seen.add(r.stage);
                return true;
              })
              .map((r) => ({
                label: r.stage,
                value: num(r.logfc) ?? 0,
                sig: (num(r.padj) ?? 1) < 0.05,
              }));
            return (
              <DivergingBar
                data={data}
                xLabel="log2 fold change"
                valueFormat={(v) => fmtLogFC(v)}
                showValues
                caption={`Axis: ${axis}${ref ? ` (vs ${ref})` : ""}`}
                height={Math.max(140, 40 + data.length * 26)}
                ariaLabel={`Stage trajectory for ${symbol}`}
              />
            );
          }}
        </QueryBlock>
      </Panel>

      {/* Sex-biased panel only renders for the canonical sex-dimorphic genes. */}
      {sex.data && sex.data.length > 0 && (
        <Panel
          title="Sex-biased expression"
          description="Female vs male effect (control-centered); this gene is in the canonical sex-dimorphic set."
        >
          {(() => {
            const r = sex.data[0];
            const data: DivergingBarRow[] = [
              { label: "Female", value: num(r.logfc_female) ?? 0 },
              { label: "Male", value: num(r.logfc_male) ?? 0 },
            ];
            return (
              <>
                <DivergingBar
                  data={data}
                  xLabel="log2 fold change"
                  valueFormat={(v) => fmtLogFC(v)}
                  showValues
                  height={120}
                  ariaLabel={`Sex-biased expression for ${symbol}`}
                />
                <p className="mt-2 text-xs text-muted-foreground">
                  Classification: {r.sex_class} · interaction padj{" "}
                  <span className="font-numeric">{fmtP(r.interaction_padj)}</span>
                </p>
              </>
            );
          })()}
        </Panel>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Genetics: colocalization across GWAS × ancestry
// ---------------------------------------------------------------------------

interface ColocRow {
  gwas: string;
  ancestry: string;
  method: string;
  pp4: number;
  pp4_susie: number | null;
  ld_reliability: string | null;
}

export function GeneticsSection({ symbol }: { symbol: string }) {
  const coloc = useParquetRows<ColocRow>(
    "coloc_by_gwas.parquet",
    `SELECT gwas, ancestry, method, pp4, pp4_susie, ld_reliability FROM coloc_by_gwas ${where(
      "symbol",
      symbol
    )} ORDER BY pp4 DESC`
  );

  const columns: DataTableColumn<ColocRow>[] = [
    { key: "gwas", header: "GWAS", sortable: true },
    {
      key: "ancestry",
      header: "Ancestry",
      render: (r) => <AncestryChip ancestry={r.ancestry} />,
      sortable: true,
    },
    { key: "method", header: "Method", sortable: true },
    {
      key: "pp4",
      header: "PP.H4",
      numeric: true,
      sortable: true,
      render: (r) => <span className="font-numeric">{fmtPP4(r.pp4)}</span>,
    },
    {
      key: "pp4_susie",
      header: "PP.H4 (SuSiE)",
      numeric: true,
      sortable: true,
      render: (r) => <span className="font-numeric">{fmtPP4(r.pp4_susie)}</span>,
    },
    { key: "ld_reliability", header: "LD reliability", sortable: true },
  ];

  return (
    <Panel
      title="Genetic colocalization"
      description="Posterior probability of a shared causal variant (PP.H4) between MASLD-related GWAS and liver eQTLs, per ancestry."
    >
      <QueryBlock q={coloc} empty="No colocalization tested for this gene.">
        {(rows) => {
          // Best PP.H4 per (gwas, ancestry) for the heatmap.
          const best = new Map<string, number>();
          for (const r of rows) {
            const k = `${r.gwas} ${r.ancestry}`;
            best.set(k, Math.max(best.get(k) ?? 0, num(r.pp4) ?? 0));
          }
          const cells: HeatmapCell[] = [...best.entries()].map(([k, v]) => {
            const [g, a] = k.split(" ");
            return { row: g, col: a, value: v };
          });
          return (
            <div className="space-y-4">
              <Heatmap
                data={cells}
                colorScale="sequential"
                domain={[0, 1]}
                showValues
                valueFormat={(v) => v.toFixed(2)}
                caption="Best PP.H4 per GWAS × ancestry (0 → 1)."
                height={Math.max(140, 60 + new Set(cells.map((c) => c.row)).size * 22)}
                ariaLabel={`Colocalization heatmap for ${symbol}`}
              />
              <DataTable
                data={rows}
                columns={columns}
                pageSize={10}
                initialSort={{ key: "pp4", dir: "desc" }}
                dense
              />
            </div>
          );
        }}
      </QueryBlock>
    </Panel>
  );
}

// ---------------------------------------------------------------------------
// Single-cell: pseudobulk DE across cell types × contrasts
// ---------------------------------------------------------------------------

interface PseudobulkRow {
  cell_type: string;
  contrast: string;
  logfc: number;
  padj: number;
}

export function SingleCellSection({ symbol }: { symbol: string }) {
  const pb = useParquetRows<PseudobulkRow>(
    "gene_pseudobulk_de.parquet",
    `SELECT cell_type, contrast, logfc, padj FROM gene_pseudobulk_de ${where(
      "symbol",
      symbol
    )}`
  );

  return (
    <Panel
      title="Cell-type resolved expression"
      description="Pseudobulk differential expression per cell type. Dot size = significance (−log10 FDR); color = log2 fold change."
    >
      <QueryBlock q={pb} empty="No single-cell pseudobulk estimates for this gene.">
        {(rows) => {
          const mag = Math.max(
            0.5,
            ...rows.map((r) => Math.abs(num(r.logfc) ?? 0))
          );
          const data: DotPlotPoint[] = rows.map((r) => ({
            row: r.cell_type,
            col: r.contrast,
            size: clamp01(-Math.log10(Math.max(num(r.padj) ?? 1, 1e-300)) / 5),
            color: num(r.logfc) ?? 0,
          }));
          return (
            <DotPlot
              data={data}
              colorDomain={[-mag, mag]}
              colorLabel="log2 FC"
              sizeLabel="−log10 FDR"
              valueFormat={(v) => fmtLogFC(v)}
              height={Math.max(160, 60 + new Set(data.map((d) => d.row)).size * 24)}
              ariaLabel={`Cell-type expression dotplot for ${symbol}`}
            />
          );
        }}
      </QueryBlock>
    </Panel>
  );
}

// ---------------------------------------------------------------------------
// Spatial: hepatic zonation
// ---------------------------------------------------------------------------

interface ZonationRow {
  zonation_class: string | null;
  spearman_rho: number | null;
  kruskal_pval: number | null;
  n_donors: number | null;
  mean_PP1: number | null;
  mean_PP2: number | null;
  mean_Mid: number | null;
  mean_PC2: number | null;
  mean_PC1: number | null;
  morans_i_healthy: number | null;
  morans_i_masld: number | null;
  delta_i: number | null;
  svg_category: string | null;
}

const ZONE_ORDER: { key: keyof ZonationRow; label: string }[] = [
  { key: "mean_PP1", label: "Periportal" },
  { key: "mean_PP2", label: "PP2" },
  { key: "mean_Mid", label: "Mid" },
  { key: "mean_PC2", label: "PC2" },
  { key: "mean_PC1", label: "Pericentral" },
];

export function SpatialSection({ symbol }: { symbol: string }) {
  const zon = useParquetRows<ZonationRow>(
    "spatial_zonation.parquet",
    `SELECT * FROM spatial_zonation ${where("symbol", symbol)}`
  );

  return (
    <Panel
      title="Hepatic zonation"
      description="Zonal expression profile from periportal (PP) to pericentral (PC) and spatial autocorrelation."
    >
      <QueryBlock q={zon} empty="This gene was not spatially profiled.">
        {(rows) => {
          const r = rows[0];
          const cells: HeatmapCell[] = ZONE_ORDER.map((z) => ({
            row: symbol,
            col: z.label,
            value: num(r[z.key]) ?? 0,
          }));
          return (
            <div className="space-y-4">
              <div className="flex flex-wrap items-center gap-2">
                {r.zonation_class && <StageChip stage={r.zonation_class} />}
                {r.svg_category && (
                  <ProgramChip label={`SVG: ${r.svg_category}`} />
                )}
              </div>
              <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
                <StatTile label="Spearman ρ" value={fmtLogFC(r.spearman_rho)} />
                <StatTile label="Moran's I (healthy)" value={fmtPP4(r.morans_i_healthy)} />
                <StatTile label="Moran's I (MASLD)" value={fmtPP4(r.morans_i_masld)} />
                <StatTile label="ΔI" value={fmtLogFC(r.delta_i)} />
              </div>
              <Heatmap
                data={cells}
                colorScale="sequential"
                showValues
                valueFormat={(v) => v.toFixed(2)}
                caption="Mean zonal expression (periportal → pericentral)."
                height={110}
                ariaLabel={`Zonation profile for ${symbol}`}
              />
            </div>
          );
        }}
      </QueryBlock>
    </Panel>
  );
}

// ---------------------------------------------------------------------------
// Programs: molecular-program membership
// ---------------------------------------------------------------------------

interface ProgramRow {
  program_code: string;
  program_label: string;
  logfc: number;
  padj: number;
  direction: string;
}

export function ProgramsSection({ symbol }: { symbol: string }) {
  const prog = useParquetRows<ProgramRow>(
    "program_gene_membership.parquet",
    `SELECT program_code, program_label, logfc, padj, direction FROM program_gene_membership ${where(
      "symbol",
      symbol
    )} ORDER BY program_code`
  );

  return (
    <Panel
      title="Molecular-program membership"
      description="Loading of this gene onto each continuous molecular-program axis (k=6)."
    >
      <QueryBlock q={prog} empty="No program membership for this gene.">
        {(rows) => {
          const data: DivergingBarRow[] = rows.map((r) => ({
            label: r.program_label || r.program_code,
            value: num(r.logfc) ?? 0,
            sig: (num(r.padj) ?? 1) < 0.05,
          }));
          return (
            <DivergingBar
              data={data}
              xLabel="program loading (log2 FC)"
              valueFormat={(v) => fmtLogFC(v)}
              showValues
              height={Math.max(140, 40 + data.length * 26)}
              ariaLabel={`Program loadings for ${symbol}`}
            />
          );
        }}
      </QueryBlock>
    </Panel>
  );
}

// ---------------------------------------------------------------------------
// Therapeutics: known drugs + reversal (LINCS) compounds
// ---------------------------------------------------------------------------

interface DrugRow {
  drug: string;
  stage: string;
  moa: string | null;
  support: string | null;
  target_class: string | null;
}
interface LincsRow {
  compound: string;
  score: number;
  moa: string | null;
}

export function TherapeuticsSection({ symbol }: { symbol: string }) {
  const drugs = useParquetRows<DrugRow>(
    "gene_drugs.parquet",
    `SELECT drug, stage, moa, support, target_class FROM gene_drugs ${where(
      "symbol",
      symbol
    )}`
  );
  const lincs = useParquetRows<LincsRow>(
    "gene_lincs.parquet",
    `SELECT compound, score, moa FROM gene_lincs ${where("symbol", symbol)} ORDER BY score`
  );

  const drugCols: DataTableColumn<DrugRow>[] = [
    { key: "drug", header: "Drug", sortable: true },
    {
      key: "stage",
      header: "Stage",
      render: (r) => <ProgramChip label={r.stage} color={devStageColor(r.stage)} />,
      sortable: true,
    },
    { key: "moa", header: "Mechanism", sortable: true },
    { key: "target_class", header: "Target class", sortable: true },
    { key: "support", header: "Support", sortable: true },
  ];
  const lincsCols: DataTableColumn<LincsRow>[] = [
    { key: "compound", header: "Compound", sortable: true },
    {
      key: "score",
      header: "Reversal score",
      numeric: true,
      sortable: true,
      render: (r) => <span className="font-numeric">{fmtLogFC(r.score)}</span>,
    },
    { key: "moa", header: "Mechanism", sortable: true },
  ];

  const noDrugs = !drugs.loading && (drugs.data?.length ?? 0) === 0;
  const noLincs = !lincs.loading && (lincs.data?.length ?? 0) === 0;
  if (noDrugs && noLincs && !drugs.error && !lincs.error) {
    return (
      <Panel title="Therapeutic actionability">
        <EmptyState
          title="No drug or reversal-compound evidence"
          description="This gene is not linked to a known drug or a signature-reversing compound in the current release."
        />
      </Panel>
    );
  }

  return (
    <div className="grid gap-4 lg:grid-cols-2">
      <Panel
        title="Drugs & clinical development"
        description="Approved or investigational agents linked to this target."
      >
        <QueryBlock q={drugs} empty="No linked drugs.">
          {(rows) => <DataTable data={rows} columns={drugCols} pageSize={10} dense />}
        </QueryBlock>
      </Panel>
      <Panel
        title="Signature-reversing compounds"
        description="LINCS compounds whose transcriptional signature opposes the disease signature (more negative = stronger reversal)."
      >
        <QueryBlock q={lincs} empty="No reversal compounds.">
          {(rows) => <DataTable data={rows} columns={lincsCols} pageSize={10} dense />}
        </QueryBlock>
      </Panel>
    </div>
  );
}
