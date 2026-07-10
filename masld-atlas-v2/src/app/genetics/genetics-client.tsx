"use client";

/**
 * /genetics — GWAS colocalization + ancestry explorer (absorbs /causal).
 *
 * Four data-backed tabs:
 *   1. COLOC Evidence  — gene × GWAS PP.H4 heatmap (coloc_by_gwas.parquet)
 *   2. Effector Genes  — SuSiE / union effector table (atlas_core.parquet)
 *   3. Ancestry        — per-gene PP.H4 forest across ancestries + portfolio
 *                        breadth + cross-ancestry-replicated genes + Broadaway
 *   4. Regulatory      — GWAS-ATAC variant browser (gwas_atac_browser.json)
 *
 * Doctrine: cross-ANCESTRY is first-class; disease claims are HUMAN-ONLY;
 * control/non-significant is `#9E9E9E`. The Regulatory tab ships a fine-mapping
 * PIP Manhattan + regulatory locus-zoom over the shipped fine-mapped variants
 * (PIP + colocalization PP.H4 + motif disruption); a genome-wide per-SNP
 * Manhattan would upgrade the same schema-forward components when that export
 * lands.
 */

import { useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";
import { PageHeader } from "@/components/page-header";
import { StatTile } from "@/components/stat-tile";
import { TabBar, type TabItem } from "@/components/tab-bar";
import { DataTable, type DataTableColumn } from "@/components/data-table";
import { AncestryChip } from "@/components/chips";
import { EmptyState, SkeletonBlock, Legend } from "@/components/states";
import {
  Heatmap,
  type HeatmapCell,
  ForestPlot,
  type ForestRow,
} from "@/components/charts";
import { queryParquet, sqlString } from "@/lib/duck";
import { dataUrl } from "@/lib/data-base";
import {
  GWAS_COUNT,
  COLOC_SUSIE,
  COLOC_UNION,
  ATLAS_GENES,
  fmt,
} from "@/lib/atlas-constants";
import { ANCESTRY_ORDER, ancestryColor } from "@/lib/palette";
import type { ManhattanVariant } from "@/components/charts/manhattan";

// Heavy per-variant viz — lazy + client-only (kept out of the shared bundle).
// Imported by DIRECT path (not the charts barrel) per the track file-ownership.
const Manhattan = dynamic(
  () => import("@/components/charts/manhattan").then((m) => m.Manhattan),
  { ssr: false, loading: () => <SkeletonBlock className="h-[320px]" /> }
);
const LocusZoom = dynamic(
  () => import("@/components/charts/locus-zoom").then((m) => m.LocusZoom),
  { ssr: false, loading: () => <SkeletonBlock className="h-[300px]" /> }
);

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------
interface ByAncestry {
  ancestry: string;
  n_gwas: number;
  n_genes_h4_05: number;
  n_genes_h4_08: number;
}
interface CrossAncGene {
  gene: string;
  best_pp4: number | null;
  ancestries: string[];
}
interface Summary {
  n_effector_genes_susie: number;
  n_effector_genes_union: number;
  n_gwas: number;
  ancestries: string[];
  by_ancestry: ByAncestry[];
  conf_tier_counts: { tier: string; n: number }[];
  cross_ancestry_replicated: CrossAncGene[];
  broadaway_overlap: { n_in_broadaway_747: number; n_novel_to_us: number };
}

interface MotifDisrupted {
  tf: string;
  effect: string;
  in_regulon: boolean;
  alleleDiff?: number;
}
interface Variant {
  id: string;
  chr: string;
  pos: number;
  pip: number;
  rec_pip?: number;
  cell_types?: string[];
  nearest_gene: string;
  distance?: number;
  motifs_disrupted?: MotifDisrupted[];
  coloc_pp4?: number;
  coloc_gwas?: string;
}
interface CellTypeEnrichment {
  cell_type: string;
  count: number;
  enrichment: number;
  padj: number;
}
interface GwasAtacData {
  variants: Variant[];
  summary: {
    total_variants: number;
    in_peaks: number;
    with_motif_disruption: number;
    disease_regulon_tfs: number;
  };
  cell_type_enrichment: CellTypeEnrichment[];
}

interface HeatRow {
  symbol: string;
  gwas: string;
  pp4: number;
}
interface EffRow {
  human_symbol: string;
  coloc_best_susie_pp4: number | null;
  coloc_best_susie_gwas: string | null;
  coloc_abf_best_pp4: number | null;
  n_coloc_sources: number | null;
  n_ancestry_gwas: number | null;
  coloc_cross_ancestry_replicated: unknown;
  coloc_susie_conf_tier: string | null;
}
interface AncPp4Row {
  ancestry: string;
  pp4: number;
}

type EffMode = "susie" | "union";
type Thr = 0.5 | 0.8;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
const DASH = "—";
const n2 = (v: number | null | undefined) =>
  v == null || Number.isNaN(v) ? DASH : v.toFixed(2);
const n3 = (v: number | null | undefined) =>
  v == null || Number.isNaN(v) ? DASH : v.toFixed(3);

function isTrue(v: unknown): boolean {
  if (typeof v === "boolean") return v;
  if (typeof v === "number") return v !== 0;
  const s = String(v).trim().toLowerCase();
  return s === "true" || s === "1" || s === "yes";
}

/** Short display label for a GWAS id, e.g. 2023_36280732_NAFLD_deCode_EUR → "NAFLD deCode (EUR)". */
function prettyGwas(g: string | null | undefined): string {
  if (!g) return DASH;
  let s = g.replace(/^\d{4}_\d+_/, ""); // strip leading YEAR_PMID_
  const anc = s.match(/_(EUR|AFR|EAS|AMR|SAS|CSA)$/);
  if (anc) s = s.slice(0, anc.index);
  s = s.replace(/_/g, " ").trim();
  return anc ? `${s} (${anc[1]})` : s;
}

const confRank = (t: string | null | undefined): number => {
  const s = (t ?? "").toLowerCase();
  return s === "high" ? 3 : s === "suggestive" ? 2 : s === "nominal" ? 1 : 0;
};
function ConfBadge({ tier }: { tier: string | null | undefined }) {
  const t = (tier ?? "").toLowerCase();
  if (t === "" || t === "none" || t === "nan")
    return <span className="text-muted-foreground">{DASH}</span>;
  const label = t.charAt(0).toUpperCase() + t.slice(1);
  const cls =
    t === "high"
      ? "border-primary/40 bg-primary/10 text-primary"
      : t === "suggestive"
        ? "border-border bg-muted text-foreground"
        : "border-border bg-muted/50 text-muted-foreground";
  return (
    <span
      className={`inline-flex rounded-full border px-2 py-0.5 text-xs font-medium ${cls}`}
    >
      {label}
    </span>
  );
}

const TABS: TabItem[] = [
  { id: "coloc", label: "COLOC Evidence" },
  { id: "effector", label: "Effector Genes" },
  { id: "ancestry", label: "Ancestry" },
  { id: "regulatory", label: "Regulatory" },
];

// ---------------------------------------------------------------------------
// Heatmap assembly (top genes × top GWAS at the chosen threshold)
// ---------------------------------------------------------------------------
const HEAT_MAX_GENES = 25;
const HEAT_MAX_GWAS = 14;

function buildHeat(rows: HeatRow[]): {
  cells: HeatmapCell[];
  rows: string[];
  cols: string[];
} {
  if (rows.length === 0) return { cells: [], rows: [], cols: [] };
  const geneBest = new Map<string, number>();
  for (const r of rows)
    geneBest.set(r.symbol, Math.max(geneBest.get(r.symbol) ?? 0, r.pp4));
  const topGenes = [...geneBest.entries()]
    .sort((a, b) => b[1] - a[1])
    .slice(0, HEAT_MAX_GENES)
    .map((e) => e[0]);
  const topGeneSet = new Set(topGenes);

  const gwasStat = new Map<string, { n: number; best: number }>();
  for (const r of rows) {
    if (!topGeneSet.has(r.symbol)) continue;
    const cur = gwasStat.get(r.gwas) ?? { n: 0, best: 0 };
    cur.n += 1;
    cur.best = Math.max(cur.best, r.pp4);
    gwasStat.set(r.gwas, cur);
  }
  const topGwas = [...gwasStat.entries()]
    .sort((a, b) => b[1].n - a[1].n || b[1].best - a[1].best)
    .slice(0, HEAT_MAX_GWAS)
    .map((e) => e[0]);
  const gwasSet = new Set(topGwas);

  // pretty labels (disambiguate the rare collision by padding a space)
  const labelOf = new Map<string, string>();
  const used = new Set<string>();
  for (const g of topGwas) {
    let lab = prettyGwas(g);
    while (used.has(lab)) lab += " ";
    used.add(lab);
    labelOf.set(g, lab);
  }

  const cells: HeatmapCell[] = [];
  for (const r of rows) {
    if (!topGeneSet.has(r.symbol) || !gwasSet.has(r.gwas)) continue;
    cells.push({ row: r.symbol, col: labelOf.get(r.gwas)!, value: r.pp4 });
  }
  return { cells, rows: topGenes, cols: topGwas.map((g) => labelOf.get(g)!) };
}

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------
export function GeneticsClient() {
  const [tab, setTab] = useState<string>("coloc");

  const [summary, setSummary] = useState<Summary | null>(null);
  const [reg, setReg] = useState<GwasAtacData | null>(null);

  const [heatThr, setHeatThr] = useState<Thr>(0.5);
  const [heatRows, setHeatRows] = useState<HeatRow[] | null>(null);

  const [effMode, setEffMode] = useState<EffMode>("susie");
  const [effRows, setEffRows] = useState<EffRow[] | null>(null);

  const [selGene, setSelGene] = useState<string>("");
  const [ancRows, setAncRows] = useState<AncPp4Row[] | null>(null);

  // ---- static JSON (summary + regulatory) --------------------------------
  useEffect(() => {
    fetch(dataUrl("coloc_ancestry_summary.json"))
      .then((r) => r.json())
      .then((d: Summary) => {
        setSummary(d);
        if (d.cross_ancestry_replicated?.length)
          setSelGene(d.cross_ancestry_replicated[0].gene);
      })
      .catch(() => setSummary(null));
    fetch(dataUrl("gwas_atac_browser.json"))
      .then((r) => r.json())
      .then((d: GwasAtacData) => setReg(d))
      .catch(() => setReg(null));
  }, []);

  // ---- COLOC heatmap (re-query on threshold change; stale-while-revalidate) ----
  useEffect(() => {
    let alive = true;
    queryParquet<HeatRow>(
      "coloc_by_gwas.parquet",
      (t) =>
        `SELECT symbol, gwas, MAX(GREATEST(pp4, COALESCE(pp4_susie, 0))) AS pp4
         FROM ${t}
         WHERE symbol <> '' AND GREATEST(pp4, COALESCE(pp4_susie, 0)) >= ${heatThr}
         GROUP BY symbol, gwas`
    )
      .then((rows) => alive && setHeatRows(rows))
      .catch(() => alive && setHeatRows([]));
    return () => {
      alive = false;
    };
  }, [heatThr]);

  // ---- effector-gene table (re-query on mode change; stale-while-revalidate) ----
  useEffect(() => {
    let alive = true;
    // Canonical Tier-1/2 effector flags (pre-computed on gene_level_coloc_tier12.csv,
    // intersected with the atlas universe by be-kickoff). Atlas-present ceilings are
    // 449 SuSiE / 975 union (the 473/1,031 headline also counts coloc-only genes
    // outside the 27,187-gene atlas — surfaced as the hero tiles + table caption).
    const where = effMode === "union" ? "coloc_tier12_union" : "coloc_tier12_pass";
    queryParquet<EffRow>(
      "atlas_core.parquet",
      (t) =>
        `SELECT human_symbol, coloc_best_susie_pp4, coloc_best_susie_gwas,
                coloc_abf_best_pp4, n_coloc_sources, n_ancestry_gwas,
                coloc_cross_ancestry_replicated, coloc_susie_conf_tier
         FROM ${t}
         WHERE ${where}
         ORDER BY GREATEST(COALESCE(coloc_best_susie_pp4, 0),
                           COALESCE(coloc_abf_best_pp4, 0)) DESC`
    )
      .then((rows) => alive && setEffRows(rows))
      .catch(() => alive && setEffRows([]));
    return () => {
      alive = false;
    };
  }, [effMode]);

  // ---- per-gene ancestry forest (re-query on gene change) ----------------
  useEffect(() => {
    if (!selGene) return;
    let alive = true;
    queryParquet<AncPp4Row>(
      "coloc_by_gwas.parquet",
      (t) =>
        `SELECT ancestry, MAX(GREATEST(pp4, COALESCE(pp4_susie, 0))) AS pp4
         FROM ${t} WHERE symbol = ${sqlString(selGene)}
         GROUP BY ancestry`
    )
      .then((rows) => alive && setAncRows(rows))
      .catch(() => alive && setAncRows([]));
    return () => {
      alive = false;
    };
  }, [selGene]);

  const heat = useMemo(() => buildHeat(heatRows ?? []), [heatRows]);

  const forest: ForestRow[] = useMemo(() => {
    if (!ancRows) return [];
    return ANCESTRY_ORDER.map((a): ForestRow | null => {
      const hit = ancRows.find((r) => r.ancestry === a);
      return hit
        ? { label: a, estimate: hit.pp4, ancestry: a, sig: hit.pp4 > 0.5 }
        : null;
    }).filter((x): x is ForestRow => x !== null);
  }, [ancRows]);

  return (
    <div className="mx-auto max-w-6xl px-6 py-8">
      <PageHeader
        eyebrow="Genetics & ancestry"
        title="Genetics & Ancestry"
        description={`Gene-level colocalization evidence and effector genes across ${GWAS_COUNT} GWAS and 5 ancestries — ${fmt(
          COLOC_SUSIE
        )} SuSiE / ${fmt(
          COLOC_UNION
        )} union effector genes. Disease claims are human-only; cross-ancestry comparisons are surfaced explicitly.`}
      />

      {/* Hero stat tiles */}
      <div className="mb-8 grid grid-cols-2 gap-3 lg:grid-cols-4">
        <StatTile label="GWAS" value={fmt(GWAS_COUNT)} sublabel="5 ancestries" />
        <StatTile label="Ancestries" value="5" sublabel="EUR · AFR · EAS · AMR · SAS" />
        <StatTile
          label="SuSiE effector genes"
          value={fmt(COLOC_SUSIE)}
          sublabel="PP.H4 > 0.5 (Tier-1/2)"
        />
        <StatTile
          label="Union effector genes"
          value={fmt(COLOC_UNION)}
          sublabel="SuSiE ∪ ABF (Tier-1/2)"
        />
      </div>

      <TabBar
        tabs={TABS}
        active={tab}
        onChange={setTab}
        aria-label="Genetics views"
        className="mb-6"
      />

      {tab === "coloc" && (
        <ColocTab heat={heat} loading={heatRows === null} thr={heatThr} setThr={setHeatThr} />
      )}
      {tab === "effector" && (
        <EffectorTab
          rows={effRows}
          mode={effMode}
          setMode={setEffMode}
        />
      )}
      {tab === "ancestry" && (
        <AncestryTab
          summary={summary}
          selGene={selGene}
          setSelGene={setSelGene}
          forest={forest}
          forestLoading={ancRows === null}
        />
      )}
      {tab === "regulatory" && <RegulatoryTab data={reg} />}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Tab 1 — COLOC Evidence
// ---------------------------------------------------------------------------
function ThresholdToggle({ thr, setThr }: { thr: Thr; setThr: (t: Thr) => void }) {
  return (
    <div className="inline-flex items-center gap-1 rounded-md border border-border p-0.5 text-xs">
      {([0.5, 0.8] as Thr[]).map((t) => (
        <button
          key={t}
          type="button"
          onClick={() => setThr(t)}
          className={`rounded px-2.5 py-1 font-medium tabular-nums transition-colors ${
            thr === t
              ? "bg-primary text-primary-foreground"
              : "text-muted-foreground hover:text-foreground"
          }`}
        >
          PP.H4 ≥ {t}
        </button>
      ))}
    </div>
  );
}

function ColocTab({
  heat,
  loading,
  thr,
  setThr,
}: {
  heat: { cells: HeatmapCell[]; rows: string[]; cols: string[] };
  loading: boolean;
  thr: Thr;
  setThr: (t: Thr) => void;
}) {
  return (
    <section>
      <div className="mb-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <h2 className="text-lg font-semibold tracking-tight">
            Gene × GWAS colocalization
          </h2>
          <p className="mt-1 max-w-2xl text-sm text-muted-foreground">
            Cell = best colocalization PP.H4 (SuSiE or ABF) per gene × GWAS. Top{" "}
            {HEAT_MAX_GENES} genes × top GWAS at the selected threshold.
          </p>
        </div>
        <ThresholdToggle thr={thr} setThr={setThr} />
      </div>

      {loading ? (
        <div className="space-y-2">
          {Array.from({ length: 12 }).map((_, i) => (
            <SkeletonBlock key={i} className="h-6" />
          ))}
        </div>
      ) : heat.cells.length === 0 ? (
        <EmptyState
          title="No colocalizations at this threshold"
          description="Lower the PP.H4 threshold to populate the heatmap."
        />
      ) : (
        <div className="overflow-x-auto rounded-lg border border-border bg-card p-4">
          <Heatmap
            data={heat.cells}
            rows={heat.rows}
            cols={heat.cols}
            colorScale="sequential"
            domain={[thr, 1]}
            valueFormat={(v) => v.toFixed(2)}
            showRowLabels
            showColLabels
            height={Math.max(340, heat.rows.length * 22 + 160)}
            ariaLabel="Gene by GWAS colocalization PP.H4 heatmap"
            caption="Gene-resolution colocalization. For variant-level fine-mapping (PIP, motif disruption), see the Regulatory tab's Manhattan + locus-zoom."
          />
        </div>
      )}
    </section>
  );
}

// ---------------------------------------------------------------------------
// Tab 2 — Effector Genes
// ---------------------------------------------------------------------------
function EffectorTab({
  rows,
  mode,
  setMode,
}: {
  rows: EffRow[] | null;
  mode: EffMode;
  setMode: (m: EffMode) => void;
}) {
  const columns: DataTableColumn<EffRow>[] = [
    { key: "human_symbol", header: "Gene", sortable: true },
    {
      key: "coloc_best_susie_pp4",
      header: "SuSiE PP.H4",
      numeric: true,
      sortable: true,
      render: (r) => n3(r.coloc_best_susie_pp4),
    },
    {
      key: "coloc_abf_best_pp4",
      header: "ABF PP.H4",
      numeric: true,
      sortable: true,
      render: (r) => n3(r.coloc_abf_best_pp4),
    },
    {
      key: "coloc_best_susie_gwas",
      header: "Best GWAS",
      sortable: true,
      sortAccessor: (r) => r.coloc_best_susie_gwas ?? "",
      render: (r) => (
        <span className="text-xs text-muted-foreground">
          {prettyGwas(r.coloc_best_susie_gwas)}
        </span>
      ),
    },
    { key: "n_coloc_sources", header: "Sources", numeric: true, sortable: true },
    { key: "n_ancestry_gwas", header: "Ancestries", numeric: true, sortable: true },
    {
      key: "coloc_cross_ancestry_replicated",
      header: "Cross-anc",
      align: "center",
      sortable: true,
      sortAccessor: (r) => (isTrue(r.coloc_cross_ancestry_replicated) ? 1 : 0),
      render: (r) =>
        isTrue(r.coloc_cross_ancestry_replicated) ? (
          <span className="text-sm font-semibold text-primary">Yes</span>
        ) : (
          <span className="text-muted-foreground">{DASH}</span>
        ),
    },
    {
      key: "coloc_susie_conf_tier",
      header: "Confidence",
      sortable: true,
      sortAccessor: (r) => confRank(r.coloc_susie_conf_tier),
      render: (r) => <ConfBadge tier={r.coloc_susie_conf_tier} />,
    },
  ];

  return (
    <section>
      <div className="mb-4 flex flex-wrap items-center justify-between gap-3">
        <div>
          <h2 className="text-lg font-semibold tracking-tight">Effector genes</h2>
          <p className="mt-1 max-w-2xl text-sm text-muted-foreground">
            Genes with a colocalizing eQTL. Toggle the SuSiE-primary set (PP.H4 &gt;
            0.5) vs the SuSiE ∪ ABF union.
          </p>
        </div>
        <div className="inline-flex items-center gap-1 rounded-md border border-border p-0.5 text-xs">
          {(["susie", "union"] as EffMode[]).map((m) => (
            <button
              key={m}
              type="button"
              onClick={() => setMode(m)}
              className={`rounded px-2.5 py-1 font-medium transition-colors ${
                mode === m
                  ? "bg-primary text-primary-foreground"
                  : "text-muted-foreground hover:text-foreground"
              }`}
            >
              {m === "susie" ? "SuSiE" : "Union"}
            </button>
          ))}
        </div>
      </div>

      {rows === null ? (
        <div className="space-y-2">
          {Array.from({ length: 10 }).map((_, i) => (
            <SkeletonBlock key={i} className="h-8" />
          ))}
        </div>
      ) : rows.length === 0 ? (
        <EmptyState title="Effector-gene table unavailable" description="atlas_core.parquet did not load." />
      ) : (
        <>
          <p className="mb-3 text-xs text-muted-foreground">
            {fmt(rows.length)} of the{" "}
            {fmt(mode === "union" ? COLOC_UNION : COLOC_SUSIE)} canonical Tier-1/2{" "}
            {mode === "union" ? "union" : "SuSiE"} effector genes appear in this
            atlas ({(mode === "union" ? COLOC_UNION : COLOC_SUSIE) - rows.length} are
            COLOC-only signals not otherwise profiled in the {fmt(ATLAS_GENES)}-gene
            atlas).
          </p>
          <DataTable
            data={rows}
            columns={columns}
            geneColumn="human_symbol"
            pageSize={20}
            initialSort={{ key: "coloc_best_susie_pp4", dir: "desc" }}
            rowKey={(r) => r.human_symbol}
          />
        </>
      )}
    </section>
  );
}

// ---------------------------------------------------------------------------
// Tab 3 — Ancestry
// ---------------------------------------------------------------------------
function AncestryTab({
  summary,
  selGene,
  setSelGene,
  forest,
  forestLoading,
}: {
  summary: Summary | null;
  selGene: string;
  setSelGene: (g: string) => void;
  forest: ForestRow[];
  forestLoading: boolean;
}) {
  if (!summary)
    return (
      <div className="space-y-2">
        {Array.from({ length: 8 }).map((_, i) => (
          <SkeletonBlock key={i} className="h-8" />
        ))}
      </div>
    );

  const maxGenes = Math.max(1, ...summary.by_ancestry.map((a) => a.n_genes_h4_05));

  const crossData = summary.cross_ancestry_replicated.map((g) => ({
    gene: g.gene,
    best_pp4: g.best_pp4,
    ancestries: g.ancestries,
  }));
  type CrossRow = (typeof crossData)[number];
  const crossCols: DataTableColumn<CrossRow>[] = [
    { key: "gene", header: "Gene", sortable: true },
    {
      key: "best_pp4",
      header: "Best PP.H4",
      numeric: true,
      sortable: true,
      render: (r) => n3(r.best_pp4),
    },
    {
      key: "ancestries",
      header: "Ancestries",
      sortable: true,
      sortAccessor: (r) => r.ancestries.length,
      render: (r) => (
        <div className="flex flex-wrap gap-1">
          {r.ancestries.map((a) => (
            <AncestryChip key={a} ancestry={a} />
          ))}
        </div>
      ),
    },
  ];

  return (
    <section className="space-y-8">
      {/* Per-gene ancestry forest */}
      <div>
        <div className="mb-3 flex flex-wrap items-center justify-between gap-3">
          <div>
            <h2 className="text-lg font-semibold tracking-tight">
              Colocalization across ancestries
            </h2>
            <p className="mt-1 max-w-2xl text-sm text-muted-foreground">
              Best PP.H4 per ancestry for a cross-ancestry-replicated gene.
              Reference line at 0.5; below-threshold points are hollow / gray.
            </p>
          </div>
          <label className="flex items-center gap-2 text-sm">
            <span className="text-muted-foreground">Gene</span>
            <select
              value={selGene}
              onChange={(e) => setSelGene(e.target.value)}
              className="rounded-md border border-border bg-background px-2 py-1 text-sm font-mono"
            >
              {summary.cross_ancestry_replicated.map((g) => (
                <option key={g.gene} value={g.gene}>
                  {g.gene}
                </option>
              ))}
            </select>
          </label>
        </div>
        <div className="rounded-lg border border-border bg-card p-4">
          {forestLoading ? (
            <SkeletonBlock className="h-40" />
          ) : forest.length === 0 ? (
            <EmptyState title="No per-ancestry colocalization for this gene" />
          ) : (
            <ForestPlot
              data={forest}
              xLabel="Best PP.H4 per ancestry"
              nullValue={0.5}
              valueFormat={(v) => v.toFixed(2)}
              height={Math.max(160, forest.length * 40 + 60)}
              ariaLabel={`PP.H4 across ancestries for ${selGene}`}
            />
          )}
        </div>
      </div>

      {/* Portfolio breadth by ancestry */}
      <div>
        <h3 className="mb-1 text-base font-semibold tracking-tight">
          Portfolio breadth by ancestry
        </h3>
        <p className="mb-3 text-sm text-muted-foreground">
          GWAS count and colocalizing genes (PP.H4 &gt; 0.5 / 0.8) per ancestry
          across the full portfolio.
        </p>
        <div className="overflow-hidden rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50 text-xs text-muted-foreground">
                <th className="px-3 py-2 text-left font-medium">Ancestry</th>
                <th className="px-3 py-2 text-right font-medium">GWAS</th>
                <th className="px-3 py-2 text-left font-medium">Genes PP.H4 &gt; 0.5</th>
                <th className="px-3 py-2 text-right font-medium">&gt; 0.8</th>
              </tr>
            </thead>
            <tbody>
              {ANCESTRY_ORDER.map((a) => {
                const row = summary.by_ancestry.find((x) => x.ancestry === a);
                if (!row) return null;
                const pct = (row.n_genes_h4_05 / maxGenes) * 100;
                return (
                  <tr key={a} className="border-b border-border/50 last:border-0">
                    <td className="px-3 py-2">
                      <AncestryChip ancestry={a} />
                    </td>
                    <td className="px-3 py-2 text-right tabular-nums">{row.n_gwas}</td>
                    <td className="px-3 py-2">
                      <div className="flex items-center gap-2">
                        <div className="h-2.5 w-32 overflow-hidden rounded-full bg-muted">
                          <div
                            className="h-2.5 rounded-full"
                            style={{
                              width: `${Math.max(pct, 2)}%`,
                              backgroundColor: ancestryColor(a),
                            }}
                          />
                        </div>
                        <span className="tabular-nums">{fmt(row.n_genes_h4_05)}</span>
                      </div>
                    </td>
                    <td className="px-3 py-2 text-right tabular-nums">
                      {fmt(row.n_genes_h4_08)}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      </div>

      {/* Broadaway overlap */}
      <div className="grid grid-cols-2 gap-3 sm:max-w-md">
        <StatTile
          label="In Broadaway-747"
          value={fmt(summary.broadaway_overlap.n_in_broadaway_747)}
          sublabel="overlap with the published liver eQTL-COLOC set"
        />
        <StatTile
          label="Novel to us"
          value={fmt(summary.broadaway_overlap.n_novel_to_us)}
          sublabel="colocalizations not in Broadaway-747"
        />
      </div>

      {/* Cross-ancestry-replicated genes */}
      <div>
        <h3 className="mb-1 text-base font-semibold tracking-tight">
          Cross-ancestry-replicated genes
        </h3>
        <p className="mb-3 text-sm text-muted-foreground">
          {fmt(summary.cross_ancestry_replicated.length)} genes colocalize in ≥ 2
          ancestries. Cross-ancestry replication reflects shared regulatory
          architecture; disease effects remain human-only.
        </p>
        <Legend
          className="mb-3"
          items={ANCESTRY_ORDER.map((a) => ({ color: ancestryColor(a), label: a }))}
        />
        <DataTable
          data={crossData}
          columns={crossCols}
          geneColumn="gene"
          pageSize={20}
          initialSort={{ key: "best_pp4", dir: "desc" }}
          rowKey={(r) => r.gene}
        />
      </div>
    </section>
  );
}

// ---------------------------------------------------------------------------
// Tab 4 — Regulatory (GWAS-ATAC)
// ---------------------------------------------------------------------------
function RegulatoryTab({ data }: { data: GwasAtacData | null }) {
  const [selectedId, setSelectedId] = useState<string | null>(null);

  // Map the fine-mapped variants onto the schema-forward Manhattan/locus-zoom
  // variant shape (motif count folded to `motifCount`).
  const manhattanVariants: ManhattanVariant[] = useMemo(
    () =>
      (data?.variants ?? []).map((v) => ({
        id: v.id,
        chr: v.chr,
        pos: v.pos,
        pip: v.pip,
        coloc_pp4: v.coloc_pp4 ?? null,
        nearest_gene: v.nearest_gene,
        motifCount: v.motifs_disrupted?.length ?? 0,
      })),
    [data]
  );

  // Default the locus-zoom to the highest-PIP variant; a Manhattan click overrides.
  const selected = useMemo(() => {
    if (manhattanVariants.length === 0) return null;
    if (selectedId) return manhattanVariants.find((v) => v.id === selectedId) ?? null;
    return [...manhattanVariants].sort((a, b) => (b.pip ?? 0) - (a.pip ?? 0))[0];
  }, [manhattanVariants, selectedId]);

  if (!data)
    return (
      <div className="space-y-2">
        {Array.from({ length: 6 }).map((_, i) => (
          <SkeletonBlock key={i} className="h-8" />
        ))}
      </div>
    );

  const varCols: DataTableColumn<Variant>[] = [
    {
      key: "id",
      header: "Variant",
      render: (r) => <span className="font-mono text-xs">{r.id}</span>,
    },
    {
      key: "locus",
      header: "Locus",
      sortable: true,
      sortAccessor: (r) => r.pos,
      render: (r) => (
        <span className="font-mono text-xs">
          {r.chr}:{r.pos.toLocaleString()}
        </span>
      ),
    },
    {
      key: "pip",
      header: "PIP",
      numeric: true,
      sortable: true,
      render: (r) => n3(r.pip),
    },
    { key: "nearest_gene", header: "Nearest gene" },
    {
      key: "coloc_pp4",
      header: "COLOC PP.H4",
      numeric: true,
      sortable: true,
      sortAccessor: (r) => r.coloc_pp4 ?? -1,
      render: (r) => n3(r.coloc_pp4),
    },
    {
      key: "motifs",
      header: "Motifs disrupted",
      numeric: true,
      sortable: true,
      sortAccessor: (r) => r.motifs_disrupted?.length ?? 0,
      render: (r) => String(r.motifs_disrupted?.length ?? 0),
    },
  ];

  return (
    <section className="space-y-8">
      <div>
        <h2 className="text-lg font-semibold tracking-tight">
          GWAS-ATAC regulatory variants
        </h2>
        <p className="mt-1 max-w-2xl text-sm text-muted-foreground">
          Fine-mapped MASLD variants intersected with liver open chromatin and
          scanned for TF-motif disruption.
        </p>
      </div>

      <div className="grid grid-cols-2 gap-3 lg:grid-cols-4">
        <StatTile label="Fine-mapped variants" value={fmt(data.summary.total_variants)} />
        <StatTile label="In ATAC peaks" value={fmt(data.summary.in_peaks)} />
        <StatTile
          label="Motif-disrupting"
          value={fmt(data.summary.with_motif_disruption)}
        />
        <StatTile
          label="Disease-regulon TFs"
          value={fmt(data.summary.disease_regulon_tfs)}
        />
      </div>

      {/* Fine-mapping PIP Manhattan + brush-linked regulatory locus-zoom */}
      {manhattanVariants.length > 0 && (
        <div className="space-y-4">
          <div>
            <h3 className="mb-1 text-base font-semibold tracking-tight">
              Fine-mapping landscape
            </h3>
            <p className="text-sm text-muted-foreground">
              Each point is a fine-mapped MASLD variant: height = posterior
              inclusion probability (PIP), color = colocalization PP.H4, and a
              regulatory-hue ring (also a tick in the rug below) marks variants
              that disrupt a TF motif. Click a point to zoom its locus.
            </p>
          </div>
          <div className="rounded-lg border border-border bg-card p-4">
            <Manhattan
              variants={manhattanVariants}
              selectedId={selected?.id ?? null}
              onVariantClick={(v) => setSelectedId(v.id)}
              height={320}
              ariaLabel="Fine-mapping PIP Manhattan of MASLD variants"
              caption="Colored by colocalization PP.H4 (gray = untested); larger points are the top-PIP lead per chromosome."
            />
          </div>
          <div className="rounded-lg border border-border bg-card p-4">
            <LocusZoom
              variants={manhattanVariants}
              selected={selected}
              onVariantClick={(v) => setSelectedId(v.id)}
              height={300}
              ariaLabel="Regulatory locus-zoom for the selected variant"
              caption={
                selected?.nearest_gene
                  ? `Locus around ${selected.nearest_gene} — nearest-gene track below the axis.`
                  : "Nearest-gene track shown below the position axis."
              }
            />
          </div>
        </div>
      )}

      {data.cell_type_enrichment?.length > 0 && (
        <div>
          <h3 className="mb-2 text-base font-semibold tracking-tight">
            Cell-type enrichment
          </h3>
          <div className="overflow-hidden rounded-lg border border-border">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border bg-muted/50 text-xs text-muted-foreground">
                  <th className="px-3 py-2 text-left font-medium">Cell type</th>
                  <th className="px-3 py-2 text-right font-medium">Variants</th>
                  <th className="px-3 py-2 text-right font-medium">Enrichment</th>
                  <th className="px-3 py-2 text-right font-medium">padj</th>
                </tr>
              </thead>
              <tbody>
                {[...data.cell_type_enrichment]
                  .sort((a, b) => b.enrichment - a.enrichment)
                  .map((c) => (
                    <tr key={c.cell_type} className="border-b border-border/50 last:border-0">
                      <td className="px-3 py-2">{c.cell_type}</td>
                      <td className="px-3 py-2 text-right tabular-nums">{c.count}</td>
                      <td className="px-3 py-2 text-right tabular-nums">
                        {n2(c.enrichment)}×
                      </td>
                      <td className="px-3 py-2 text-right tabular-nums">
                        {c.padj < 0.001 ? c.padj.toExponential(1) : n3(c.padj)}
                      </td>
                    </tr>
                  ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      <div>
        <h3 className="mb-2 text-base font-semibold tracking-tight">
          Top fine-mapped variants
        </h3>
        <DataTable
          data={data.variants}
          columns={varCols}
          geneColumn="nearest_gene"
          pageSize={20}
          initialSort={{ key: "pip", dir: "desc" }}
          rowKey={(r, i) => `${r.id}-${i}`}
        />
      </div>
    </section>
  );
}
