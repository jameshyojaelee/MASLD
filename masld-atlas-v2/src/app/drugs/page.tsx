"use client";

import { useEffect, useMemo, useState } from "react";
import { PageHeader } from "@/components/page-header";
import { PageContainer } from "@/components/page-container";
import { StatTile } from "@/components/stat-tile";
import { DataTable, type DataTableColumn } from "@/components/data-table";
import { TabBar, type TabItem } from "@/components/tab-bar";
import {
  StackedBar,
  Heatmap,
  Bar,
  type HeatmapCell,
  type BarDatum,
} from "@/components/charts";
import { EmptyState, SkeletonBlock, Legend } from "@/components/states";
import { Badge } from "@/components/ui/badge";
import { queryParquet } from "@/lib/duck";
import { dataUrl } from "@/lib/data-base";
import { DEV_STAGE_COLORS, CONTROL, MODALITY_HEX } from "@/lib/palette";
import {
  DRUGS_APPROVED,
  DRUGS_CLINICAL,
  DRUGS_PRECLINICAL,
  fmt,
} from "@/lib/atlas-constants";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface TopCompound {
  rank: number;
  name: string;
  target: string | null;
  moa: string | null;
  composite_score: number;
  score_reversal: number;
  score_network: number;
}

interface GeneDrug {
  symbol: string;
  drug: string;
  stage: string;
  moa: string;
  support: string;
}

interface GeneLincs {
  symbol: string;
  compound: string;
  score: number;
  moa: string;
}

interface PhaseCount {
  stage: string;
  count: number;
}

interface HeatRow {
  human_symbol: string;
  bulk_shrunk_logFC: number | null;
  coloc_best_susie_pp4: number | null;
  coloc_abf_best_pp4: number | null;
  twas_z: number | null;
  essentiality_chronos: number | null;
  spatial_morans_i: number | null;
  pharos_tdl: string | null;
  max_phase_masld: number | null;
}

// Supplementary essentiality section: most essential *druggable* atlas genes by
// DepMap CHRONOS gene-effect. `essentiality_chronos` is guaranteed non-null by
// the query's WHERE clause (so a plain `number` here is safe).
interface EssentialGene {
  human_symbol: string;
  essentiality_chronos: number;
  pharos_tdl: string | null;
  max_phase_masld: number | null;
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

// Phase → color. Approved/Clinical from the sanctioned dev-stage palette;
// Preclinical uses a light-blue extension of the same ramp rather than the
// control gray (#9E9E9E is reserved for control/healthy series, never a data
// series like "preclinical").
const PHASE_ORDER = ["Approved", "Clinical", "Preclinical"] as const;
const PHASE_COLOR: Record<string, string> = {
  Approved: DEV_STAGE_COLORS.Approved,
  Clinical: DEV_STAGE_COLORS.Clinical,
  Preclinical: "#a9cce3",
};

// Map the drug-target drug_dev_status buckets onto the three headline
// development phases. Non-MASLD buckets (undetermined / discovery /
// drugged_other_indication / discontinued) are the broader druggable universe,
// not pipeline phases, so they are excluded from the phase decomposition.
// Source = drug_targets.parquet (full drug-target universe, NOT atlas-restricted)
// so the GROUP BY lands the canonical 2 / 208 / 1,504 exactly.
const DEV_STATUS_TO_PHASE: Record<string, string> = {
  masld_approved: "Approved",
  masld_clinical: "Clinical",
  masld_preclinical: "Preclinical",
};

// Fallback phase counts (drug-level headline from atlas-constants), used when
// the live parquet query is unavailable.
const FALLBACK_PHASES: PhaseCount[] = [
  { stage: "Approved", count: DRUGS_APPROVED },
  { stage: "Clinical", count: DRUGS_CLINICAL },
  { stage: "Preclinical", count: DRUGS_PRECLINICAL },
];

// Evidence-modality columns for the target × evidence heatmap. Each returns a
// value normalised to [0, 1] (or null → blank cell) so modalities on different
// native scales are visually comparable.
const clamp01 = (x: number) => (x < 0 ? 0 : x > 1 ? 1 : x);

function pharosScore(tdl: string | null): number | null {
  switch (tdl) {
    case "Tclin":
      return 1;
    case "Tchem":
      return 0.7;
    case "Tbio":
      return 0.4;
    case "Tdark":
      return 0.15;
    default:
      return null;
  }
}

const EVIDENCE_COLUMNS: {
  col: string;
  get: (r: HeatRow) => number | null;
}[] = [
  {
    col: "Human DE",
    get: (r) =>
      r.bulk_shrunk_logFC == null ? null : clamp01(Math.abs(r.bulk_shrunk_logFC) / 2),
  },
  {
    col: "COLOC",
    get: (r) => {
      const v = Math.max(r.coloc_best_susie_pp4 ?? 0, r.coloc_abf_best_pp4 ?? 0);
      return v > 0 ? clamp01(v) : null;
    },
  },
  {
    col: "TWAS",
    get: (r) => (r.twas_z == null ? null : clamp01(Math.abs(r.twas_z) / 6)),
  },
  {
    col: "Essentiality",
    get: (r) =>
      r.essentiality_chronos == null ? null : clamp01(-r.essentiality_chronos),
  },
  {
    col: "Spatial",
    get: (r) =>
      r.spatial_morans_i == null ? null : clamp01(r.spatial_morans_i / 0.5),
  },
  {
    col: "Druggability",
    get: (r) => {
      const phase = r.max_phase_masld ? r.max_phase_masld / 4 : null;
      const pharos = pharosScore(r.pharos_tdl);
      const best = Math.max(phase ?? 0, pharos ?? 0);
      return best > 0 ? clamp01(best) : null;
    },
  },
];

// ---------------------------------------------------------------------------
// Page
// ---------------------------------------------------------------------------

export default function DrugsPage() {
  const [phases, setPhases] = useState<PhaseCount[] | null>(null);
  const [druggableCount, setDruggableCount] = useState<number | null>(null);
  const [heatRows, setHeatRows] = useState<HeatRow[]>([]);
  const [essentialGenes, setEssentialGenes] = useState<EssentialGene[]>([]);
  const [topCompounds, setTopCompounds] = useState<TopCompound[]>([]);
  const [geneDrugs, setGeneDrugs] = useState<GeneDrug[]>([]);
  const [geneLincs, setGeneLincs] = useState<GeneLincs[]>([]);
  const [loading, setLoading] = useState(true);
  const [tableTab, setTableTab] = useState("reversal");

  useEffect(() => {
    let cancelled = false;

    async function load() {
      // Phase decomposition (live from the full drug-target universe, so the
      // canonical 2 / 208 / 1,504 is reproduced exactly — not atlas-restricted).
      const phasePromise = queryParquet<{ status: string; n: number }>(
        "drug_targets.parquet",
        (t) =>
          `SELECT drug_dev_status AS status, COUNT(*)::INT AS n
           FROM ${t} WHERE drug_dev_status IS NOT NULL
           GROUP BY drug_dev_status`
      )
        .then((rows) => {
          const acc: Record<string, number> = {
            Approved: 0,
            Clinical: 0,
            Preclinical: 0,
          };
          for (const r of rows) {
            const phase = DEV_STATUS_TO_PHASE[r.status];
            if (phase) acc[phase] += r.n;
          }
          const out = PHASE_ORDER.map((s) => ({ stage: s, count: acc[s] }));
          return out.some((p) => p.count > 0) ? out : FALLBACK_PHASES;
        })
        .catch(() => FALLBACK_PHASES);

      const druggablePromise = queryParquet<{ n: number }>(
        "atlas_core.parquet",
        (t) => `SELECT COUNT(*)::INT AS n FROM ${t} WHERE dgidb_druggable`
      )
        .then((rows) => rows[0]?.n ?? null)
        .catch(() => null);

      // Target × evidence heatmap rows: DE-supported druggable targets, ranked
      // by breadth of evidence (layers_active), then convergence.
      const heatPromise = queryParquet<HeatRow>(
        "atlas_core.parquet",
        (t) =>
          `SELECT human_symbol, bulk_shrunk_logFC,
                  coloc_best_susie_pp4, coloc_abf_best_pp4, twas_z,
                  essentiality_chronos, spatial_morans_i, pharos_tdl, max_phase_masld
           FROM ${t}
           WHERE is_deg
             AND (dgidb_druggable
                  OR pharos_tdl IN ('Tclin','Tchem')
                  OR max_phase_masld > 0
                  OR drug_dev_status LIKE 'masld_%')
           ORDER BY layers_active DESC, convergence_score DESC,
                    ABS(bulk_shrunk_logFC) DESC
           LIMIT 24`
      ).catch(() => [] as HeatRow[]);

      // Supplementary: most essential *druggable* atlas genes by DepMap CHRONOS
      // gene-effect (more negative = more essential). Not a convergence channel —
      // essentiality is a constitutive property, reported here as context only.
      const essentialPromise = queryParquet<EssentialGene>(
        "atlas_core.parquet",
        (t) =>
          `SELECT human_symbol, essentiality_chronos, pharos_tdl, max_phase_masld
           FROM ${t}
           WHERE dgidb_druggable AND essentiality_chronos IS NOT NULL
           ORDER BY essentiality_chronos ASC
           LIMIT 15`
      ).catch(() => [] as EssentialGene[]);

      // Reversal compounds (JSON, refreshed by the data pipeline).
      const compoundsPromise = fetch(dataUrl("drug_pipeline.json"))
        .then((r) => (r.ok ? r.json() : null))
        .then((d) => {
          const raw: unknown[] = Array.isArray(d?.top_compounds)
            ? d.top_compounds
            : [];
          return raw.map((c, i) => {
            const o = c as Record<string, unknown>;
            return {
              rank: i + 1,
              name: String(o.name ?? ""),
              target: (o.target as string | null) ?? null,
              moa: (o.moa as string | null) ?? null,
              composite_score: Number(o.composite_score ?? 0),
              score_reversal: Number(o.score_reversal ?? 0),
              score_network: Number(o.score_network ?? 0),
            } satisfies TopCompound;
          });
        })
        .catch(() => [] as TopCompound[]);

      // Per-target clinical drugs + LINCS reversal compounds.
      const drugsPromise = queryParquet<GeneDrug>(
        "gene_drugs.parquet",
        (t) => `SELECT symbol, drug, stage, moa, support FROM ${t}`
      ).catch(() => [] as GeneDrug[]);

      const lincsPromise = queryParquet<GeneLincs>(
        "gene_lincs.parquet",
        (t) => `SELECT symbol, compound, score, moa FROM ${t} ORDER BY score DESC`
      ).catch(() => [] as GeneLincs[]);

      const [ph, dg, heat, ess, comp, drugs, lincs] = await Promise.all([
        phasePromise,
        druggablePromise,
        heatPromise,
        essentialPromise,
        compoundsPromise,
        drugsPromise,
        lincsPromise,
      ]);

      if (cancelled) return;
      setPhases(ph);
      setDruggableCount(dg);
      setHeatRows(heat);
      setEssentialGenes(ess);
      setTopCompounds(comp);
      setGeneDrugs(drugs);
      setGeneLincs(lincs);
      setLoading(false);
    }

    load();
    return () => {
      cancelled = true;
    };
  }, []);

  // ----- derived: phase bar (long format) --------------------------------
  const phaseData = useMemo(
    () =>
      (phases ?? []).map((p) => ({
        group: p.stage,
        category: p.stage,
        value: p.count,
      })),
    [phases]
  );

  const phaseCount = (stage: string) =>
    phases?.find((p) => p.stage === stage)?.count ?? 0;

  // ----- derived: heatmap cells ------------------------------------------
  const heatData = useMemo<HeatmapCell[]>(() => {
    const cells: HeatmapCell[] = [];
    for (const r of heatRows) {
      for (const { col, get } of EVIDENCE_COLUMNS) {
        const v = get(r);
        if (v != null && !Number.isNaN(v)) {
          cells.push({ row: r.human_symbol, col, value: v });
        }
      }
    }
    return cells;
  }, [heatRows]);

  const heatRowOrder = useMemo(() => heatRows.map((r) => r.human_symbol), [heatRows]);

  // ----- derived: essentiality bars (magnitude = −CHRONOS) ----------------
  // Data mark colored with the sanctioned essentiality-modality hex from the
  // palette (never a raw utility class). More-negative CHRONOS → longer bar.
  const essentialBars = useMemo<BarDatum[]>(
    () =>
      essentialGenes.map((g) => ({
        label: g.human_symbol,
        value: -g.essentiality_chronos,
        color: MODALITY_HEX.s3_essential,
        annotation: g.essentiality_chronos.toFixed(2),
      })),
    [essentialGenes]
  );

  // ----- table column defs -----------------------------------------------
  const compoundCols: DataTableColumn<TopCompound>[] = [
    { key: "rank", header: "#", numeric: true, width: "3rem" },
    { key: "name", header: "Compound", sortable: true },
    { key: "target", header: "Target" },
    { key: "moa", header: "Mechanism" },
    {
      key: "composite_score",
      header: "Composite",
      numeric: true,
      sortable: true,
      render: (r) => r.composite_score.toFixed(3),
    },
    {
      key: "score_reversal",
      header: "Reversal",
      numeric: true,
      sortable: true,
      render: (r) => r.score_reversal.toFixed(3),
    },
    {
      key: "score_network",
      header: "Network",
      numeric: true,
      sortable: true,
      render: (r) => (r.score_network > 0 ? r.score_network.toFixed(3) : "—"),
    },
  ];

  const geneDrugCols: DataTableColumn<GeneDrug>[] = [
    { key: "symbol", header: "Target", sortable: true },
    { key: "drug", header: "Drug", sortable: true },
    { key: "stage", header: "Stage", sortable: true },
    { key: "moa", header: "Mechanism" },
    { key: "support", header: "Atlas support", sortable: true },
  ];

  const lincsCols: DataTableColumn<GeneLincs>[] = [
    { key: "symbol", header: "Target", sortable: true },
    { key: "compound", header: "Compound", sortable: true },
    { key: "moa", header: "Mechanism" },
    {
      key: "score",
      header: "Reversal score",
      numeric: true,
      sortable: true,
      render: (r) => r.score.toFixed(3),
    },
  ];

  const essentialCols: DataTableColumn<EssentialGene>[] = [
    { key: "human_symbol", header: "Gene", sortable: true },
    {
      key: "essentiality_chronos",
      header: "CHRONOS",
      numeric: true,
      sortable: true,
      render: (r) => r.essentiality_chronos.toFixed(3),
    },
    {
      key: "pharos_tdl",
      header: "Druggable · TDL",
      sortable: true,
      render: (r) => (r.pharos_tdl ? `DGIdb · ${r.pharos_tdl}` : "DGIdb"),
    },
    {
      key: "max_phase_masld",
      header: "MASLD max phase",
      numeric: true,
      sortable: true,
      render: (r) =>
        r.max_phase_masld != null && r.max_phase_masld > 0
          ? String(r.max_phase_masld)
          : "—",
    },
  ];

  const tableTabs: TabItem[] = [
    { id: "reversal", label: `Top reversal compounds (${topCompounds.length})` },
    { id: "clinical", label: `Clinical MASLD drugs (${geneDrugs.length})` },
    { id: "lincs", label: `Reversal by target (${geneLincs.length})` },
  ];

  return (
    <PageContainer>
      <PageHeader
        eyebrow="Therapeutics"
        title="Drug Pipeline"
        description={
          <>
            Translating human MASLD transcriptomic signatures into therapeutic
            hypotheses through LINCS L1000 signature reversal, network proximity,
            druggability, and clinical-stage validation. The MASLD drug landscape
            spans {DRUGS_APPROVED} approved, {fmt(DRUGS_CLINICAL)} clinical-stage,
            and {fmt(DRUGS_PRECLINICAL)} preclinical programs.
          </>
        }
      />

      {/* KPI tiles */}
      <div className="mb-10 grid grid-cols-2 gap-3 lg:grid-cols-4">
        <StatTile label="Approved" value={fmt(DRUGS_APPROVED)} sublabel="for MASLD/MASH" />
        <StatTile
          label="In clinical trials"
          value={fmt(DRUGS_CLINICAL)}
          sublabel="clinical-stage programs"
        />
        <StatTile
          label="Preclinical"
          value={fmt(DRUGS_PRECLINICAL)}
          sublabel="discovery / preclinical"
        />
        <StatTile
          label="Druggable genes"
          value={druggableCount != null ? fmt(druggableCount) : "—"}
          sublabel="DGIdb druggable genome"
        />
      </div>

      {/* Development-phase decomposition (retires the funnel) */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          MASLD drug-development stage
        </h2>
        <p className="mb-4 max-w-3xl text-sm text-muted-foreground">
          MASLD drug targets grouped by the most-advanced development stage of a
          drug acting on them. The pipeline is bottom-heavy: a large preclinical
          space narrows to a handful of clinical-stage and approved therapies.
        </p>
        <div className="rounded-lg border border-border bg-card p-4 shadow-sm">
          {loading ? (
            <SkeletonBlock className="h-[300px]" />
          ) : phaseData.length > 0 ? (
            <>
              <StackedBar
                data={phaseData}
                groups={[...PHASE_ORDER]}
                categories={[...PHASE_ORDER]}
                colorFor={(c) => PHASE_COLOR[c] ?? CONTROL}
                yLabel="targets"
                showLegend={false}
                height={300}
                ariaLabel="MASLD drug targets by development stage"
              />
              <Legend
                className="mt-2 justify-center"
                items={PHASE_ORDER.map((s) => ({
                  color: PHASE_COLOR[s],
                  label: `${s} (${fmt(phaseCount(s))})`,
                }))}
              />
            </>
          ) : (
            <EmptyState title="No development-stage data available." />
          )}
        </div>
      </section>

      {/* Target × evidence heatmap */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Target × evidence
        </h2>
        <p className="mb-4 max-w-3xl text-sm text-muted-foreground">
          Druggable disease genes scored across independent evidence modalities.
          Each column is normalised to a comparable 0–1 scale; blank cells mark
          modalities with no measurement for that target.
        </p>
        <div className="rounded-lg border border-border bg-card p-4 shadow-sm">
          {loading ? (
            <SkeletonBlock className="h-[420px]" />
          ) : heatData.length > 0 ? (
            <Heatmap
              data={heatData}
              rows={heatRowOrder}
              cols={EVIDENCE_COLUMNS.map((c) => c.col)}
              colorScale="sequential"
              domain={[0, 1]}
              valueFormat={(v) => v.toFixed(2)}
              ariaLabel="Druggable target by evidence-modality heatmap"
              height={Math.max(320, heatRowOrder.length * 22 + 90)}
            />
          ) : (
            <EmptyState title="No druggable-target evidence available." />
          )}
        </div>
      </section>

      {/* Reversal / clinical / per-target tables */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Compounds &amp; drugs
        </h2>
        <p className="mb-4 max-w-3xl text-sm text-muted-foreground">
          Signature-reversal compounds ranked by composite score, clinical-stage
          MASLD drugs with atlas transcriptomic support, and per-target LINCS
          reversal hits.
        </p>

        <TabBar
          tabs={tableTabs}
          active={tableTab}
          onChange={setTableTab}
          aria-label="Compound tables"
          className="mb-4"
        />

        {loading ? (
          <div className="space-y-2">
            <SkeletonBlock className="h-8" />
            <SkeletonBlock className="h-64" />
          </div>
        ) : tableTab === "reversal" ? (
          topCompounds.length > 0 ? (
            <DataTable
              data={topCompounds}
              columns={compoundCols}
              geneColumn="target"
              rowKey={(r) => `${r.rank}-${r.name}`}
              initialSort={{ key: "composite_score", dir: "desc" }}
              pageSize={25}
            />
          ) : (
            <EmptyState title="Reversal-compound data unavailable." />
          )
        ) : tableTab === "clinical" ? (
          geneDrugs.length > 0 ? (
            <DataTable
              data={geneDrugs}
              columns={geneDrugCols}
              geneColumn="symbol"
              rowKey={(r) => `${r.symbol}-${r.drug}`}
              pageSize={25}
            />
          ) : (
            <EmptyState title="Clinical-drug data unavailable." />
          )
        ) : geneLincs.length > 0 ? (
          <DataTable
            data={geneLincs}
            columns={lincsCols}
            geneColumn="symbol"
            rowKey={(r) => `${r.symbol}-${r.compound}`}
            initialSort={{ key: "score", dir: "desc" }}
            pageSize={25}
          />
        ) : (
          <EmptyState title="LINCS reversal data unavailable." />
        )}
      </section>

      {/* Essentiality (DepMap) — supplementary, NOT a convergence modality */}
      <section className="mb-10">
        <div className="mb-1 flex items-center gap-2">
          <Badge
            variant="secondary"
            className="uppercase tracking-wide"
          >
            Supplementary
          </Badge>
          <h2 className="text-xl font-semibold tracking-tight">
            Essentiality (DepMap CHRONOS)
          </h2>
        </div>
        <p className="mb-4 max-w-3xl text-sm text-muted-foreground">
          DepMap CHRONOS gene-effect scores flag genes whose knockout reduces
          cancer-cell-line viability — a druggability / target-tractability
          signal (more negative = more essential). Essentiality is a{" "}
          <em>constitutive</em> cellular property, not MASLD-specific evidence,
          so it is reported here as supplementary context rather than a
          convergence channel. Shown: the most essential druggable atlas genes
          (DGIdb-druggable).
        </p>
        <div className="grid gap-4 lg:grid-cols-2">
          <div className="rounded-lg border border-border bg-card p-4 shadow-sm">
            {loading ? (
              <SkeletonBlock className="h-[380px]" />
            ) : essentialBars.length > 0 ? (
              <Bar
                data={essentialBars}
                orientation="horizontal"
                valueLabel="−CHRONOS (more essential →)"
                height={Math.max(320, essentialBars.length * 22 + 60)}
                valueFormat={(v) => v.toFixed(2)}
                ariaLabel="Most essential druggable atlas genes by DepMap CHRONOS gene-effect"
                tooltipLines={(d) => (
                  <>
                    <div className="font-medium">{d.label}</div>
                    <div>CHRONOS {(-d.value).toFixed(3)}</div>
                  </>
                )}
              />
            ) : (
              <EmptyState title="No essentiality data available." />
            )}
          </div>
          <div className="rounded-lg border border-border bg-card p-4 shadow-sm">
            {loading ? (
              <SkeletonBlock className="h-[380px]" />
            ) : essentialGenes.length > 0 ? (
              <DataTable
                data={essentialGenes}
                columns={essentialCols}
                geneColumn="human_symbol"
                rowKey={(r) => r.human_symbol}
                initialSort={{ key: "essentiality_chronos", dir: "asc" }}
                pageSize={0}
                dense
              />
            ) : (
              <EmptyState title="No essentiality data available." />
            )}
          </div>
        </div>
      </section>
    </PageContainer>
  );
}
