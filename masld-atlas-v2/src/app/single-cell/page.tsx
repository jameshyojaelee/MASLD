"use client";

/**
 * /single-cell — multi-scale single-cell + spatial view over the 1.23M-cell
 * integrated atlas. Four tabs: Composition, UMAP, Cell-type DE, Spatial.
 *
 * Data:
 *   singlecell_summary.json      composition, per-cell-type disease activity,
 *                                fGSEA F-transitions, hepatocyte subtypes
 *   sc_umap_downsampled.parquet  ~60k-cell embedding (regl-scatterplot)
 *   gene_pseudobulk_de.parquet   per-cell-type DE (contrast-picked)
 *   sc_pseudobulk_de.parquet     all-cell pseudobulk DE
 *   sc_hep_markers.parquet       hepatocyte subtype markers (dot plot)
 *   spatial_summary.json         zonation + SVG summaries (single-cohort)
 *   spatial_zonation.parquet     per-gene zone-bin means (heatmap)
 */

import { useEffect, useMemo, useState, type ReactNode } from "react";
import dynamic from "next/dynamic";
import { PageHeader } from "@/components/page-header";
import { StatTile } from "@/components/stat-tile";
import { TabBar } from "@/components/tab-bar";
import { GeneLink } from "@/components/gene-link";
import { EmptyState, SkeletonBlock } from "@/components/states";
import { StackedBar, DivergingBar, Volcano, DotPlot, Heatmap } from "@/components/charts";
import { dataUrl } from "@/lib/data-base";
import { queryParquet } from "@/lib/duck";
import { fmt, SC_CELLS, DEG_GATE_LABEL } from "@/lib/atlas-constants";
import type { ColorMode, UmapPoint } from "./umap-scatter";

const UmapScatter = dynamic(() => import("./umap-scatter"), {
  ssr: false,
  loading: () => (
    <div className="space-y-3">
      {/* Square, shape-matched skeleton so the embedding doesn't pop in. */}
      <div className="relative h-[520px] w-full overflow-hidden rounded-lg border border-border bg-card">
        <div className="absolute inset-0 animate-pulse bg-muted/40" />
        <div className="absolute left-1/2 top-1/2 h-40 w-40 -translate-x-1/2 -translate-y-1/2 animate-pulse rounded-full bg-muted/60 blur-xl" />
      </div>
      <div className="flex flex-wrap gap-x-4 gap-y-1.5">
        {Array.from({ length: 6 }).map((_, i) => (
          <SkeletonBlock key={i} className="h-3 w-20" />
        ))}
      </div>
    </div>
  ),
});

// ---------------------------------------------------------------------------
// JSON shapes (subset of singlecell_summary.json / spatial_summary.json)
// ---------------------------------------------------------------------------

interface CompositionRow {
  sample: string;
  dataset: string;
  condition: string;
  condition_harmonized: string;
  fractions: Record<string, number | null>;
}
interface HepSubtype {
  subtype: string;
  n_cells: number;
  odds_ratio: number | null;
  fisher_padj: number | null;
  enrichment_class: string;
  stage_spearman_rho: number | null;
  axes: Record<string, number | null>;
}
interface FgseaRow {
  cell_type: string;
  contrast: string;
  nes: number | null;
  padj: number | null;
  size: number;
  leading_edge_n: number;
  significant: boolean;
}
interface ScSummary {
  n_cells: number;
  n_datasets: number;
  n_samples: number;
  cell_types: string[];
  composition: CompositionRow[];
  celltype_disease_activity: { cell_type: string; activity: number | null }[];
  fgsea_transitions: FgseaRow[];
  hep_subtypes: HepSubtype[];
}
interface SpatialSummary {
  n_samples: number;
  datasets: string[];
  zonation_dist: { zonation_bin: string; n: number }[];
  zonation_by_condition: { condition: string; mean_zonation_score: number | null; n: number }[];
  top_disrupted: {
    gene: string;
    rho_healthy: number | null;
    rho_masld: number | null;
    disruption_score: number | null;
    underpowered: boolean;
  }[];
  svg_categories: { category: string; n: number }[];
  top_differential_svgs: {
    gene: string;
    morans_i_healthy: number | null;
    morans_i_masld: number | null;
    delta_i: number | null;
    category: string;
  }[];
}

type TabId = "composition" | "umap" | "de" | "spatial";

const NBSP_DASH = "—";

// ---------------------------------------------------------------------------
// Page shell
// ---------------------------------------------------------------------------

export default function SingleCellPage() {
  const [summary, setSummary] = useState<ScSummary | null>(null);
  const [spatial, setSpatial] = useState<SpatialSummary | null>(null);
  const [tab, setTab] = useState<TabId>("composition");
  // Cross-tab jump: a UMAP lasso can pre-set the Cell-type DE tab's cell type.
  const [dePreset, setDePreset] = useState<string | null>(null);

  useEffect(() => {
    fetch(dataUrl("singlecell_summary.json"))
      .then((r) => r.json())
      .then(setSummary)
      .catch(() => setSummary(null));
    fetch(dataUrl("spatial_summary.json"))
      .then((r) => r.json())
      .then(setSpatial)
      .catch(() => setSpatial(null));
  }, []);

  return (
    <div className="mx-auto max-w-6xl px-6 py-8">
      <PageHeader
        eyebrow="Multi-scale"
        title="Single-Cell & Spatial"
        description={`Cell-type composition, a ${fmt(
          SC_CELLS
        )}-cell UMAP embedding, per-cell-type differential expression, and spatial zonation across the integrated liver atlas.`}
      />

      {/* Hero stats */}
      <div className="mb-8 grid grid-cols-2 gap-4 sm:grid-cols-4">
        <StatTile label="Cells" value={summary ? fmt(summary.n_cells) : fmt(SC_CELLS)} />
        <StatTile label="Samples" value={summary ? fmt(summary.n_samples) : NBSP_DASH} />
        <StatTile
          label="Cell types"
          value={summary ? summary.cell_types.length : NBSP_DASH}
        />
        <StatTile label="Datasets" value={summary ? summary.n_datasets : NBSP_DASH} />
      </div>

      <TabBar
        aria-label="Single-cell views"
        active={tab}
        onChange={(id) => setTab(id as TabId)}
        tabs={[
          { id: "composition", label: "Composition" },
          { id: "umap", label: "UMAP" },
          { id: "de", label: "Cell-type DE" },
          { id: "spatial", label: "Spatial" },
        ]}
        className="mb-6"
      />

      {tab === "composition" && <CompositionTab summary={summary} />}
      {tab === "umap" && (
        <UmapTab
          onViewDE={(cellType) => {
            setDePreset(cellType);
            setTab("de");
          }}
        />
      )}
      {tab === "de" && (
        <CellTypeDeTab cellTypes={summary?.cell_types ?? []} initialCellType={dePreset} />
      )}
      {tab === "spatial" && <SpatialTab spatial={spatial} />}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Tab 1 — Composition
// ---------------------------------------------------------------------------

function CompositionTab({ summary }: { summary: ScSummary | null }) {
  const [dataset, setDataset] = useState<string>("All datasets");

  const datasets = useMemo(() => {
    if (!summary) return [] as string[];
    return ["All datasets", ...Array.from(new Set(summary.composition.map((c) => c.dataset))).sort()];
  }, [summary]);

  const bars = useMemo(() => {
    if (!summary) return [];
    const rows = summary.composition.filter(
      (c) => dataset === "All datasets" || c.dataset === dataset
    );
    // Mean cell-type fraction per disease state (condition_harmonized).
    const groups = new Map<string, { sum: Record<string, number>; n: number }>();
    for (const r of rows) {
      const g = r.condition_harmonized;
      if (!groups.has(g)) groups.set(g, { sum: {}, n: 0 });
      const acc = groups.get(g)!;
      acc.n += 1;
      for (const ct of summary.cell_types) acc.sum[ct] = (acc.sum[ct] ?? 0) + (r.fractions[ct] ?? 0);
    }
    const out: { group: string; category: string; value: number }[] = [];
    for (const [g, acc] of groups) {
      for (const ct of summary.cell_types)
        out.push({ group: g, category: ct, value: acc.n ? acc.sum[ct] / acc.n : 0 });
    }
    return out;
  }, [summary, dataset]);

  const activity = useMemo(() => {
    if (!summary) return [];
    return [...summary.celltype_disease_activity]
      .filter((a) => a.activity != null)
      .sort((a, b) => (b.activity ?? 0) - (a.activity ?? 0))
      .map((a) => ({ label: a.cell_type, value: a.activity as number }));
  }, [summary]);

  if (!summary) return <LoadingPanel />;

  const groupOrder = ["Healthy", "MASLD"].filter((g) => bars.some((b) => b.group === g));

  return (
    <div className="space-y-10">
      <section>
        <div className="mb-3 flex flex-wrap items-center justify-between gap-3">
          <div>
            <h2 className="text-lg font-semibold">Cell-type composition by disease state</h2>
            <p className="mt-1 text-sm text-muted-foreground">
              Mean cell-type fraction per sample, aggregated across{" "}
              {dataset === "All datasets" ? `${summary.n_datasets} datasets` : dataset}. Bars
              normalize to 100%.
            </p>
          </div>
          <Picker value={dataset} onChange={setDataset} options={datasets} label="Dataset" />
        </div>
        <StackedBar
          data={bars}
          groups={groupOrder}
          categories={summary.cell_types}
          normalize
          yLabel="fraction of cells"
          xLabel="disease state"
          height={360}
          ariaLabel="Cell-type composition stacked bar"
        />
      </section>

      <section>
        <h2 className="text-lg font-semibold">Per-cell-type disease activity</h2>
        <p className="mb-3 mt-1 text-sm text-muted-foreground">
          Signed MASLD disease-signature activity per cell type (positive = higher in disease).
          Derived from the cell-type disease-activity matrix.
        </p>
        {activity.length > 0 ? (
          <DivergingBar
            data={activity}
            xLabel="disease-signature activity (z)"
            clamp={1.5}
            showValues
            ariaLabel="Per-cell-type disease activity"
          />
        ) : (
          <EmptyState title="No activity scores available" />
        )}
      </section>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Tab 2 — UMAP
// ---------------------------------------------------------------------------

interface UmapRow {
  cell_id: string;
  umap_x: number;
  umap_y: number;
  cell_type: string;
  dataset: string;
  condition_harmonized: string;
}

function UmapTab({ onViewDE }: { onViewDE: (cellType: string) => void }) {
  const [points, setPoints] = useState<UmapPoint[] | null>(null);
  const [error, setError] = useState(false);
  const [colorMode, setColorMode] = useState<ColorMode>("cell_type");
  const [reduce, setReduce] = useState(false);
  const [lasso, setLasso] = useState(false);
  const [selection, setSelection] = useState<UmapPoint[]>([]);
  const [clearNonce, setClearNonce] = useState(0);

  useEffect(() => {
    let alive = true;
    queryParquet<UmapRow>(
      "sc_umap_downsampled.parquet",
      (t) =>
        `SELECT cell_id, umap_x, umap_y, cell_type, dataset, condition_harmonized FROM ${t}`
    )
      .then((rows) => {
        if (!alive) return;
        setPoints(
          rows.map((r) => ({
            x: r.umap_x,
            y: r.umap_y,
            cellType: r.cell_type,
            dataset: r.dataset,
            condition: r.condition_harmonized,
            cellId: r.cell_id,
          }))
        );
      })
      .catch(() => alive && setError(true));
    return () => {
      alive = false;
    };
  }, []);

  // Client-side tallies over the lassoed cells (counts only — raw DE is not
  // computed on arbitrary selections; the DE jump routes to the pooled cell-type
  // pseudobulk contrast instead).
  const summary = useMemo(() => {
    if (selection.length === 0) return null;
    const tally = (key: (p: UmapPoint) => string) => {
      const m = new Map<string, number>();
      for (const p of selection) m.set(key(p), (m.get(key(p)) ?? 0) + 1);
      return [...m.entries()].sort((a, b) => b[1] - a[1]);
    };
    const byCellType = tally((p) => p.cellType);
    return {
      n: selection.length,
      byCellType,
      byCondition: tally((p) => p.condition),
      byDataset: tally((p) => p.dataset),
      dominant: byCellType[0]?.[0] ?? "",
    };
  }, [selection]);

  if (error) return <EmptyState title="Failed to load UMAP embedding" />;
  if (!points) return <LoadingPanel label="Loading 60k-cell embedding…" />;

  return (
    <div className="space-y-4">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <div>
          <h2 className="text-lg font-semibold">UMAP embedding</h2>
          <p className="mt-1 text-sm text-muted-foreground">
            {fmt(points.length)} cells (downsampled, stratified by cell type × disease state).
          </p>
        </div>
        <div className="flex flex-wrap items-center gap-2">
          <Segmented
            value={colorMode}
            onChange={(v) => setColorMode(v as ColorMode)}
            options={[
              { id: "cell_type", label: "Cell type" },
              { id: "condition_harmonized", label: "Disease state" },
            ]}
          />
          <Toggle active={lasso} onClick={() => setLasso((l) => !l)}>
            Lasso
          </Toggle>
          <Toggle active={reduce} onClick={() => setReduce((r) => !r)}>
            Reduce points
          </Toggle>
        </div>
      </div>

      <div className="grid gap-4 lg:grid-cols-[1fr_260px]">
        <UmapScatter
          points={points}
          colorMode={colorMode}
          reduce={reduce}
          lasso={lasso}
          onSelect={setSelection}
          clearNonce={clearNonce}
        />

        {/* Selection summary panel */}
        <aside className="rounded-lg border border-border bg-card p-4">
          {!summary ? (
            <div className="text-sm text-muted-foreground">
              <p className="font-medium text-foreground">No selection</p>
              <p className="mt-1">
                Turn on <span className="font-medium">Lasso</span> and drag a region
                to tally its cells and jump to that cell type&apos;s differential
                expression.
              </p>
            </div>
          ) : (
            <div className="space-y-3 text-sm">
              <div className="flex items-baseline justify-between">
                <span className="font-medium text-foreground">Selection</span>
                <span className="tabular-nums text-muted-foreground">
                  {fmt(summary.n)} cells
                </span>
              </div>
              <SummaryList title="Cell type" rows={summary.byCellType} />
              <SummaryList title="Disease state" rows={summary.byCondition} />
              <SummaryList title="Dataset" rows={summary.byDataset} max={4} />
              <div className="flex flex-col gap-2 pt-1">
                {summary.dominant && (
                  <button
                    type="button"
                    onClick={() => onViewDE(summary.dominant)}
                    className="rounded-md bg-primary px-3 py-1.5 text-sm font-medium text-primary-foreground transition-opacity hover:opacity-90"
                  >
                    View DE for {summary.dominant.replace(/_/g, " ")}
                  </button>
                )}
                <button
                  type="button"
                  onClick={() => {
                    setSelection([]);
                    setClearNonce((n) => n + 1);
                  }}
                  className="rounded-md border border-border px-3 py-1.5 text-sm text-muted-foreground transition-colors hover:text-foreground"
                >
                  Clear selection
                </button>
              </div>
            </div>
          )}
        </aside>
      </div>

      <p className="text-xs text-muted-foreground">
        Hover a cell for its type, disease state, and dataset. Toggle{" "}
        <span className="font-medium">Lasso</span>, then drag to marquee-select; scroll to zoom.
      </p>
    </div>
  );
}

/** Compact ranked count list for the UMAP selection summary. */
function SummaryList({
  title,
  rows,
  max = 6,
}: {
  title: string;
  rows: [string, number][];
  max?: number;
}) {
  const shown = rows.slice(0, max);
  const rest = rows.length - shown.length;
  return (
    <div>
      <div className="mb-1 text-xs font-medium uppercase tracking-wide text-muted-foreground">
        {title}
      </div>
      <ul className="space-y-0.5">
        {shown.map(([k, v]) => (
          <li key={k} className="flex items-center justify-between gap-2">
            <span className="truncate">{k.replace(/_/g, " ")}</span>
            <span className="shrink-0 tabular-nums text-muted-foreground">{fmt(v)}</span>
          </li>
        ))}
        {rest > 0 && (
          <li className="text-xs text-muted-foreground">+{rest} more</li>
        )}
      </ul>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Tab 3 — Cell-type DE
// ---------------------------------------------------------------------------

const ALL_CELLS = "All cells";
const CONTRASTS = ["disease_vs_control", "Steatosis_vs_Healthy", "Steatohepatitis_vs_Steatosis"];

interface DeRow {
  symbol: string;
  logfc: number;
  padj: number;
}
interface HepMarkerRow {
  symbol: string;
  subtype: string;
  logfoldchanges: number;
  pvals_adj: number;
  pct_nz_group: number;
}

function CellTypeDeTab({
  cellTypes,
  initialCellType,
}: {
  cellTypes: string[];
  initialCellType?: string | null;
}) {
  const [cellType, setCellType] = useState<string>(initialCellType ?? ALL_CELLS);
  const [contrast, setContrast] = useState<string>("disease_vs_control");
  const [rows, setRows] = useState<DeRow[] | null>(null);
  const [markers, setMarkers] = useState<HepMarkerRow[] | null>(null);

  // Cell-type options: All cells + the per-cell-type DE cell types.
  const [deCellTypes, setDeCellTypes] = useState<string[]>([]);
  useEffect(() => {
    queryParquet<{ cell_type: string }>(
      "gene_pseudobulk_de.parquet",
      (t) => `SELECT DISTINCT cell_type FROM ${t} ORDER BY cell_type`
    )
      .then((r) => setDeCellTypes(r.map((x) => x.cell_type)))
      .catch(() => setDeCellTypes([]));
  }, []);

  // Load the selected DE set.
  useEffect(() => {
    let alive = true;
    setRows(null);
    const load =
      cellType === ALL_CELLS
        ? queryParquet<DeRow>(
            "sc_pseudobulk_de.parquet",
            (t) =>
              `SELECT symbol, logfc, padj FROM ${t} WHERE padj IS NOT NULL AND logfc IS NOT NULL`
          )
        : queryParquet<DeRow>(
            "gene_pseudobulk_de.parquet",
            (t) =>
              `SELECT symbol, logfc, padj FROM ${t} WHERE cell_type = ? AND contrast = ? AND padj IS NOT NULL AND logfc IS NOT NULL`,
            [cellType, contrast]
          );
    load.then((r) => alive && setRows(r)).catch(() => alive && setRows([]));
    return () => {
      alive = false;
    };
  }, [cellType, contrast]);

  // Hepatocyte subtype markers (first 10 subtypes, top markers each).
  useEffect(() => {
    queryParquet<HepMarkerRow>(
      "sc_hep_markers.parquet",
      (t) => `SELECT symbol, subtype, logfoldchanges, pvals_adj, pct_nz_group FROM ${t}`
    )
      .then((r) => setMarkers(r))
      .catch(() => setMarkers([]));
  }, []);

  const volcanoData = useMemo(() => {
    if (!rows) return [];
    return rows.map((r) => ({ symbol: r.symbol, logFC: r.logfc, padj: r.padj }));
  }, [rows]);

  const highlight = useMemo(() => {
    if (!rows) return [];
    return [...rows]
      .filter((r) => r.padj < 0.05)
      .sort((a, b) => Math.abs(b.logfc) - Math.abs(a.logfc))
      .slice(0, 8)
      .map((r) => r.symbol);
  }, [rows]);

  const nSig = useMemo(() => (rows ? rows.filter((r) => r.padj < 0.05).length : 0), [rows]);

  const dotData = useMemo(() => {
    if (!markers) return { points: [] as { row: string; col: string; size: number; color: number }[], rows: [] as string[], cols: [] as string[] };
    const subs = Array.from(new Set(markers.map((m) => m.subtype)))
      .sort((a, b) => Number(a) - Number(b))
      .slice(0, 10);
    const subSet = new Set(subs);
    const perSub = new Map<string, HepMarkerRow[]>();
    for (const m of markers) {
      if (!subSet.has(m.subtype)) continue;
      const arr = perSub.get(m.subtype) ?? [];
      arr.push(m);
      perSub.set(m.subtype, arr);
    }
    const geneOrder: string[] = [];
    for (const s of subs) {
      const top = (perSub.get(s) ?? [])
        .sort((a, b) => b.logfoldchanges - a.logfoldchanges)
        .slice(0, 3);
      for (const g of top) if (!geneOrder.includes(g.symbol)) geneOrder.push(g.symbol);
    }
    const geneSet = new Set(geneOrder);
    const cols = subs.map((s) => `H${s}`);
    const points = markers
      .filter((m) => subSet.has(m.subtype) && geneSet.has(m.symbol))
      .map((m) => ({
        row: m.symbol,
        col: `H${m.subtype}`,
        size: Math.max(0, Math.min(1, m.pct_nz_group)),
        color: m.logfoldchanges,
      }));
    return { points, rows: geneOrder, cols };
  }, [markers]);

  return (
    <div className="space-y-10">
      <section>
        <div className="mb-3 flex flex-wrap items-end justify-between gap-3">
          <div>
            <h2 className="text-lg font-semibold">Pseudobulk differential expression</h2>
            <p className="mt-1 text-sm text-muted-foreground">
              {cellType === ALL_CELLS
                ? "All-cell pseudobulk DE (disease vs control)."
                : `${cellType.replace(/_/g, " ")} — ${contrast.replace(/_/g, " ")}.`}{" "}
              {rows ? `${fmt(nSig)} genes at FDR < 0.05.` : ""}
            </p>
          </div>
          <div className="flex flex-wrap items-center gap-2">
            <Picker
              value={cellType}
              onChange={setCellType}
              options={[ALL_CELLS, ...(deCellTypes.length ? deCellTypes : cellTypes)]}
              label="Cell type"
            />
            {cellType !== ALL_CELLS && (
              <Picker value={contrast} onChange={setContrast} options={CONTRASTS} label="Contrast" />
            )}
          </div>
        </div>
        {!rows ? (
          <LoadingPanel />
        ) : rows.length === 0 ? (
          <EmptyState title="No differential expression rows for this selection" />
        ) : (
          <Volcano
            data={volcanoData}
            highlightSymbols={highlight}
            enableBrush
            height={400}
            caption={DEG_GATE_LABEL}
            ariaLabel="Cell-type DE volcano"
          />
        )}
        {rows && rows.length > 0 && (
          <div className="mt-3 flex flex-wrap gap-x-4 gap-y-1 text-xs text-muted-foreground">
            <span>Top genes:</span>
            {highlight.map((s) => (
              <GeneLink key={s} symbol={s} />
            ))}
          </div>
        )}
      </section>

      <section>
        <h2 className="text-lg font-semibold">Hepatocyte subtype markers</h2>
        <p className="mb-3 mt-1 text-sm text-muted-foreground">
          Marker genes across hepatocyte subclusters (H0–H9). Dot area = fraction of cells
          expressing; color = log fold-change vs other subtypes.
        </p>
        {!markers ? (
          <LoadingPanel />
        ) : dotData.points.length === 0 ? (
          <EmptyState title="No marker data available" />
        ) : (
          <DotPlot
            data={dotData.points}
            rows={dotData.rows}
            cols={dotData.cols}
            colorLabel="log FC"
            sizeLabel="% expressing"
            ariaLabel="Hepatocyte subtype marker dot plot"
          />
        )}
      </section>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Tab 4 — Spatial
// ---------------------------------------------------------------------------

interface ZonationRow {
  symbol: string;
  zonation_class: string;
  spearman_rho: number | null;
  mean_PP1: number | null;
  mean_PP2: number | null;
  mean_Mid: number | null;
  mean_PC2: number | null;
  mean_PC1: number | null;
}

const ZONE_COLS: { key: keyof ZonationRow; label: string }[] = [
  { key: "mean_PP1", label: "PP1" },
  { key: "mean_PP2", label: "PP2" },
  { key: "mean_Mid", label: "Mid" },
  { key: "mean_PC2", label: "PC2" },
  { key: "mean_PC1", label: "PC1" },
];

function SpatialTab({ spatial }: { spatial: SpatialSummary | null }) {
  const [zon, setZon] = useState<ZonationRow[] | null>(null);

  useEffect(() => {
    queryParquet<ZonationRow>(
      "spatial_zonation.parquet",
      (t) =>
        `SELECT symbol, zonation_class, spearman_rho, mean_PP1, mean_PP2, mean_Mid, mean_PC2, mean_PC1
         FROM ${t}
         WHERE spearman_rho IS NOT NULL AND mean_PP1 IS NOT NULL
         ORDER BY ABS(spearman_rho) DESC
         LIMIT 24`
    )
      .then((r) => setZon(r))
      .catch(() => setZon([]));
  }, []);

  // Per-gene min-max scaled zone means → shows the periportal→pericentral pattern.
  const heat = useMemo(() => {
    if (!zon) return [];
    const cells: { row: string; col: string; value: number }[] = [];
    for (const g of zon) {
      const vals = ZONE_COLS.map((c) => (g[c.key] as number | null) ?? null);
      const nums = vals.filter((v): v is number => v != null);
      if (nums.length === 0) continue;
      const lo = Math.min(...nums);
      const hi = Math.max(...nums);
      const span = hi - lo || 1;
      ZONE_COLS.forEach((c, i) => {
        const v = vals[i];
        if (v != null) cells.push({ row: g.symbol, col: c.label, value: (v - lo) / span });
      });
    }
    return cells;
  }, [zon]);

  const rowOrder = useMemo(
    () =>
      zon
        ? [...zon]
            .sort((a, b) => (b.spearman_rho ?? 0) - (a.spearman_rho ?? 0))
            .map((g) => g.symbol)
        : [],
    [zon]
  );

  return (
    <div className="space-y-8">
      {/* Single-cohort caveat */}
      <div className="rounded-lg border border-amber-500/40 bg-amber-50/60 px-4 py-3 text-sm text-amber-900 dark:bg-amber-950/30 dark:text-amber-200">
        Spatial zonation is derived from a <strong>single Visium cohort</strong>
        {spatial ? ` (${spatial.datasets.join(", ")}, ${spatial.n_samples} samples)` : ""} with a
        known batch confound. Patterns are shown <strong>descriptively</strong> (periportal ↔
        pericentral organization) and are <strong>not</strong> a disease-direction claim.
      </div>

      {/* Summary tiles */}
      <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
        <StatTile label="Visium samples" value={spatial ? spatial.n_samples : NBSP_DASH} />
        <StatTile
          label="SVGs profiled"
          value={spatial ? fmt(spatial.svg_categories.reduce((a, c) => a + c.n, 0)) : NBSP_DASH}
        />
        <StatTile
          label="Emergent SVGs"
          value={
            spatial
              ? fmt(spatial.svg_categories.find((c) => c.category === "disease_emergent_SVG")?.n ?? 0)
              : NBSP_DASH
          }
          sublabel="gain spatial structure"
        />
        <StatTile
          label="Lost SVGs"
          value={
            spatial
              ? fmt(spatial.svg_categories.find((c) => c.category === "disease_lost_SVG")?.n ?? 0)
              : NBSP_DASH
          }
          sublabel="lose spatial structure"
        />
      </div>

      {/* Zonation heatmap */}
      <section>
        <h2 className="text-lg font-semibold">Zonation profiles of the most zonated genes</h2>
        <p className="mb-3 mt-1 text-sm text-muted-foreground">
          Top 24 genes by |Spearman ρ| along the porto-central axis. Each row is min–max scaled
          across zones, so color shows the periportal (PP1) → pericentral (PC1) gradient, not
          absolute expression.
        </p>
        {!zon ? (
          <LoadingPanel />
        ) : heat.length === 0 ? (
          <EmptyState title="No zonation data available" />
        ) : (
          <Heatmap
            data={heat}
            rows={rowOrder}
            cols={ZONE_COLS.map((c) => c.label)}
            colorScale="sequential"
            valueFormat={(v) => v.toFixed(2)}
            height={Math.max(220, rowOrder.length * 20 + 90)}
            ariaLabel="Gene by zone-bin zonation heatmap"
          />
        )}
      </section>

      {/* Top disrupted genes */}
      {spatial && spatial.top_disrupted.length > 0 && (
        <section>
          <h2 className="text-lg font-semibold">Genes with the most disrupted zonation</h2>
          <p className="mb-3 mt-1 text-sm text-muted-foreground">
            Largest shift in porto-central correlation between healthy and MASLD spots
            (descriptive; underpowered genes flagged).
          </p>
          <div className="overflow-x-auto rounded-lg border border-border">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border bg-muted/50 text-xs text-muted-foreground">
                  <th className="px-3 py-2 text-left font-medium">Gene</th>
                  <th className="px-3 py-2 text-right font-medium">ρ healthy</th>
                  <th className="px-3 py-2 text-right font-medium">ρ MASLD</th>
                  <th className="px-3 py-2 text-right font-medium">disruption</th>
                  <th className="px-3 py-2 text-left font-medium">flag</th>
                </tr>
              </thead>
              <tbody>
                {spatial.top_disrupted.slice(0, 12).map((r) => (
                  <tr key={r.gene} className="border-b border-border/50 hover:bg-muted/30">
                    <td className="px-3 py-1.5">
                      <GeneLink symbol={r.gene} />
                    </td>
                    <td className="px-3 py-1.5 text-right font-mono text-xs">
                      {r.rho_healthy != null ? r.rho_healthy.toFixed(2) : NBSP_DASH}
                    </td>
                    <td className="px-3 py-1.5 text-right font-mono text-xs">
                      {r.rho_masld != null ? r.rho_masld.toFixed(2) : NBSP_DASH}
                    </td>
                    <td className="px-3 py-1.5 text-right font-mono text-xs">
                      {r.disruption_score != null ? r.disruption_score.toFixed(3) : NBSP_DASH}
                    </td>
                    <td className="px-3 py-1.5 text-xs text-muted-foreground">
                      {r.underpowered ? "underpowered" : ""}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Small shared UI
// ---------------------------------------------------------------------------

function LoadingPanel({ label }: { label?: string }) {
  return (
    <div className="space-y-3 rounded-lg border border-border bg-card p-6">
      {label && <p className="text-sm text-muted-foreground">{label}</p>}
      <SkeletonBlock className="h-6 w-1/3" />
      <SkeletonBlock className="h-48" />
    </div>
  );
}

function Picker({
  value,
  onChange,
  options,
  label,
}: {
  value: string;
  onChange: (v: string) => void;
  options: string[];
  label: string;
}) {
  return (
    <label className="flex items-center gap-2 text-xs text-muted-foreground">
      <span className="font-medium uppercase tracking-wide">{label}</span>
      <select
        value={value}
        onChange={(e) => onChange(e.target.value)}
        className="rounded-md border border-border bg-background px-2 py-1 text-sm text-foreground outline-none focus-visible:ring-2 focus-visible:ring-ring/50"
      >
        {options.map((o) => (
          <option key={o} value={o}>
            {o.replace(/_/g, " ")}
          </option>
        ))}
      </select>
    </label>
  );
}

function Segmented({
  value,
  onChange,
  options,
}: {
  value: string;
  onChange: (v: string) => void;
  options: { id: string; label: string }[];
}) {
  return (
    <div className="inline-flex overflow-hidden rounded-md border border-border">
      {options.map((o) => (
        <button
          key={o.id}
          type="button"
          onClick={() => onChange(o.id)}
          className={`px-3 py-1 text-sm transition-colors ${
            value === o.id
              ? "bg-primary text-primary-foreground"
              : "bg-background text-muted-foreground hover:text-foreground"
          }`}
        >
          {o.label}
        </button>
      ))}
    </div>
  );
}

function Toggle({
  active,
  onClick,
  children,
}: {
  active: boolean;
  onClick: () => void;
  children: ReactNode;
}) {
  return (
    <button
      type="button"
      onClick={onClick}
      aria-pressed={active}
      className={`rounded-md border px-3 py-1 text-sm transition-colors ${
        active
          ? "border-primary bg-primary/10 text-foreground"
          : "border-border bg-background text-muted-foreground hover:text-foreground"
      }`}
    >
      {children}
    </button>
  );
}
