"use client";

/**
 * /programs — interactive client.
 *
 * Two tabs:
 *  1. Programs   — k=6 NMF programs (P1–P6) as CONTINUOUS interpretive axes
 *                  (2026-07-04 audit: NOT discrete patient subtypes), each with
 *                  its curated top-gene loading bars + a program×stage activity
 *                  heatmap, plus the 117 Hotspot autocorrelation modules and the
 *                  cNMF change-points.
 *  2. Sex-biased — the canonical 8 sex-dimorphic genes (7 male / 1 female) under
 *                  the LVQW-fixed (C2) sex layer. "Sex-biased" framing ONLY; the
 *                  female-enriched-progressor narrative is retired (2026-06-08).
 *
 * Data: programs_summary.json + program_gene_membership.parquet (duck) +
 * sex_summary.json. All gene marks link to #/gene?symbol= via GeneLink.
 *
 * NOTE on the loading bars: raw top-|logFC| in program_gene_membership.parquet
 * is dominated by sex-chromosome escapees (KDM5D/USP9Y/XIST) that the atlas NMF
 * confounder-strip removes. We therefore restrict each program's bars to its
 * curated `top_genes` (already confounder-aware) and read their per-program
 * logFC from the parquet, ordered by |logFC|.
 */

import { useEffect, useMemo, useState } from "react";
import { Check } from "lucide-react";
import { TabBar } from "@/components/tab-bar";
import { DataTable, type DataTableColumn } from "@/components/data-table";
import { StatTile } from "@/components/stat-tile";
import { GeneLink } from "@/components/gene-link";
import { ProgramChip } from "@/components/chips";
import { SkeletonBlock, EmptyState, Legend } from "@/components/states";
import {
  DivergingBar,
  Heatmap,
  type DivergingBarRow,
  type HeatmapCell,
} from "@/components/charts";
import { dataUrl } from "@/lib/data-base";
import { queryParquet, sqlString } from "@/lib/duck";
import { sequentialColor, CONTROL, DIVERGING } from "@/lib/palette";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface StageActivity {
  stage: string;
  mean_score: number;
}
interface NmfProgram {
  code: string;
  label: string;
  label_category: string;
  top_pathway: string;
  top_genes: string[];
  mean_activity_by_stage: StageActivity[];
}
interface HotspotModule {
  cell_type: string;
  module: string;
  module_name: string;
  module_hallmark: string;
  is_novel: boolean;
  stability_score: number;
  bulk_replicated: boolean;
  n_genes: number;
  mean_bulk_logfc: number;
  progression_module: boolean;
}
interface Changepoint {
  program: string;
  best_model: string;
  breakpoint_stage: number | null;
  delta_aic: number;
}
interface ProgramsSummary {
  nmf_programs: NmfProgram[];
  hotspot_modules: HotspotModule[];
  n_hotspot_modules: number;
  n_novel: number;
  changepoints: Changepoint[];
}
interface MembershipRow {
  program_code: string;
  symbol: string;
  logfc: number;
  padj: number | null;
  direction: string;
}
interface SexGene {
  symbol: string;
  ensembl: string;
  chr: string;
  biotype: string;
  sex_class: string;
  logfc_female: number;
  logfc_male: number;
  interaction_padj: number;
  beta_female: number;
  beta_male: number;
}
interface SexSummary {
  n_dimorphic: number;
  n_female_biased: number;
  n_male_biased: number;
  genes: SexGene[];
  note: string;
}
interface SpatialProgramCoverageRow {
  release_id: string;
  dataset_id: string;
  assay_id: string;
  program_uid: string;
  program_label: string;
  n_program_genes: number;
  n_genes_measured: number;
  retained_l1_weight: number;
  coverage_status: string;
  testability_reason: string;
  dataset_gate: string;
  biological_unit: string;
  biological_unit_resolution: string;
  source_dependence: string;
}
interface SpatialProgramEffectRow {
  release_id: string;
  dataset_id: string;
  assay_id: string;
  program_uid: string;
  estimand: string;
  effect_unit: string;
  estimate: number | null;
  qvalue: number | null;
  n_biological: number | null;
  n_technical: number | null;
  evidence_state: string;
  testability_reason: string;
}

// ---------------------------------------------------------------------------
// Constants + small helpers
// ---------------------------------------------------------------------------

const STAGES = ["F0", "F1", "F2", "F3", "F4"];
const TOP_GENES_PER_PROGRAM = 10;

// Categorical program-identity colors (dots only; never overrides control gray).
const PROGRAM_COLORS: Record<string, string> = {
  P1: "#4C78A8",
  P2: "#F58518",
  P3: "#54A24B",
  P4: "#B279A2",
  P5: "#72B7B2",
  P6: "#E45756",
};
const programColor = (code: string): string => PROGRAM_COLORS[code] ?? CONTROL;

// Sex-bias identity colors (categorical dot only).
const SEX_COLORS: Record<string, string> = {
  Male_biased: "#3B6FB6",
  Female_biased: "#C4569A",
};

function prettySexClass(s: string): string {
  return s.replace(/_/g, "-");
}

/** Split a "HALLMARK_FOO_BAR(NES=2.8,padj=…)" label into a clean name + NES. */
function formatPathway(raw: string): { name: string; nes: string | null } {
  if (!raw) return { name: "—", nes: null };
  const m = raw.match(/^(.*?)\((.*)\)\s*$/);
  const rawName = (m ? m[1] : raw).trim();
  const inside = m ? m[2] : null;
  const name = rawName
    .replace(/^(HALLMARK|GOBP|GOCC|GOMF|REACTOME|KEGG|WP)_/, "")
    .replace(/_/g, " ")
    .toLowerCase()
    .replace(/\b\w/g, (c) => c.toUpperCase());
  let nes: string | null = null;
  if (inside) {
    const hit = inside.match(/NES\s*=\s*(-?[0-9.]+)/i);
    if (hit) nes = hit[1];
  }
  return { name: name || "—", nes };
}

function fmtSignedLFC(v: number): string {
  return v >= 0 ? `+${v.toFixed(2)}` : v.toFixed(2);
}

function fmtInteractionPadj(v: number): string {
  if (v === 0) return "<1e-16";
  return v.toExponential(1);
}

// ---------------------------------------------------------------------------
// Sub-components
// ---------------------------------------------------------------------------

/** A compact F0–F4 activity strip (one program), colored on the global scale. */
function StageStrip({
  stages,
  domain,
}: {
  stages: StageActivity[];
  domain: [number, number];
}) {
  const [lo, hi] = domain;
  const span = hi - lo || 1;
  return (
    <div className="flex items-end gap-1.5">
      {stages.map((s) => {
        const t = (s.mean_score - lo) / span;
        return (
          <div key={s.stage} className="flex flex-col items-center gap-1">
            <div
              className="h-6 w-9 rounded-sm border border-border/60"
              style={{ backgroundColor: sequentialColor(t) }}
              title={`${s.stage}: ${s.mean_score.toFixed(3)}`}
            />
            <span className="text-[10px] tabular-nums text-muted-foreground">
              {s.stage}
            </span>
          </div>
        );
      })}
    </div>
  );
}

function ProgramCard({
  program,
  bars,
  scoreDomain,
}: {
  program: NmfProgram;
  bars: MembershipRow[];
  scoreDomain: [number, number];
}) {
  const pathway = formatPathway(program.top_pathway);
  const barData: DivergingBarRow[] = bars.map((r) => ({
    label: r.symbol,
    value: r.logfc,
    sig: r.padj == null ? true : r.padj < 0.05,
  }));
  const clamp = Math.max(
    1,
    Math.ceil(Math.max(0, ...barData.map((b) => Math.abs(b.value))))
  );

  return (
    <div className="flex flex-col gap-3 rounded-lg border border-border bg-card p-4 shadow-sm">
      <div className="flex flex-wrap items-center justify-between gap-2">
        <div className="flex items-center gap-2">
          <ProgramChip label={program.code} color={programColor(program.code)} />
          <span className="text-sm font-semibold text-foreground">
            {program.label}
          </span>
        </div>
        <span className="text-xs text-muted-foreground">
          {program.label_category}
        </span>
      </div>

      <p className="text-xs text-muted-foreground">
        Top pathway: <span className="text-foreground">{pathway.name}</span>
        {pathway.nes != null && (
          <span className="tabular-nums"> · NES {pathway.nes}</span>
        )}
      </p>

      {barData.length > 0 ? (
        <DivergingBar
          data={barData}
          xLabel="program marker log₂FC"
          clamp={clamp}
          showValues
          ariaLabel={`${program.code} top program genes`}
        />
      ) : (
        <p className="text-xs text-muted-foreground">
          No curated program genes resolved.
        </p>
      )}

      <div className="flex items-center justify-between gap-3 border-t border-border/60 pt-3">
        <span className="text-[11px] font-medium uppercase tracking-wider text-muted-foreground">
          Activity by stage
        </span>
        <StageStrip stages={program.mean_activity_by_stage} domain={scoreDomain} />
      </div>
    </div>
  );
}

/** Boolean flag rendered as a check (present) / em-dash (absent) — icon, not color. */
function BoolFlag({ value }: { value: boolean }) {
  return value ? (
    <Check className="mx-auto size-4 text-foreground" aria-label="yes" />
  ) : (
    <span className="text-muted-foreground">—</span>
  );
}

// ---------------------------------------------------------------------------
// Programs tab
// ---------------------------------------------------------------------------

function ProgramsTab({
  data,
  topGenes,
}: {
  data: ProgramsSummary;
  topGenes: Record<string, MembershipRow[]>;
}) {
  const programs = data.nmf_programs;

  const allScores = programs.flatMap((p) =>
    p.mean_activity_by_stage.map((s) => s.mean_score)
  );
  const scoreDomain: [number, number] = [
    Math.min(...allScores),
    Math.max(...allScores),
  ];

  const heatmapCells: HeatmapCell[] = programs.flatMap((p) =>
    p.mean_activity_by_stage.map((s) => ({
      row: p.code,
      col: s.stage,
      value: s.mean_score,
    }))
  );

  const hotspotCols: DataTableColumn<HotspotModule>[] = [
    {
      key: "cell_type",
      header: "Cell type",
      sortable: true,
      render: (r) => <span className="capitalize">{r.cell_type}</span>,
    },
    { key: "module_name", header: "Module", sortable: true },
    {
      key: "module_hallmark",
      header: "Hallmark",
      render: (r) => (
        <span className="text-xs text-muted-foreground">{r.module_hallmark}</span>
      ),
    },
    { key: "n_genes", header: "Genes", numeric: true, sortable: true },
    {
      key: "stability_score",
      header: "Stability",
      numeric: true,
      sortable: true,
      render: (r) => r.stability_score.toFixed(2),
    },
    {
      key: "is_novel",
      header: "Novel",
      align: "center",
      sortable: true,
      sortAccessor: (r) => (r.is_novel ? 1 : 0),
      render: (r) =>
        r.is_novel ? (
          <span className="inline-flex rounded-full border border-border bg-muted px-2 py-0.5 text-[11px] font-medium">
            Novel
          </span>
        ) : (
          <span className="text-muted-foreground">—</span>
        ),
    },
    {
      key: "bulk_replicated",
      header: "Bulk-replicated",
      align: "center",
      sortable: true,
      sortAccessor: (r) => (r.bulk_replicated ? 1 : 0),
      render: (r) => <BoolFlag value={r.bulk_replicated} />,
    },
    {
      key: "progression_module",
      header: "Progression",
      align: "center",
      sortable: true,
      sortAccessor: (r) => (r.progression_module ? 1 : 0),
      render: (r) => <BoolFlag value={r.progression_module} />,
    },
  ];

  const changepointCols: DataTableColumn<Changepoint>[] = [
    {
      key: "program",
      header: "cNMF program",
      sortable: true,
      render: (r) => (
        <span className="font-mono text-xs">
          {r.program.replace(/^cnmf_global_k16_/, "k16 · ")}
        </span>
      ),
    },
    { key: "best_model", header: "Best model", sortable: true },
    {
      key: "breakpoint_stage",
      header: "Breakpoint",
      numeric: true,
      sortable: true,
      render: (r) =>
        r.breakpoint_stage == null ? (
          <span className="text-muted-foreground">—</span>
        ) : (
          r.breakpoint_stage.toFixed(2)
        ),
    },
    {
      key: "delta_aic",
      header: "ΔAIC",
      numeric: true,
      sortable: true,
      render: (r) => r.delta_aic.toFixed(2),
    },
  ];

  return (
    <div className="space-y-10">
      {/* Summary tiles */}
      <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
        <StatTile label="NMF programs" value={programs.length} sublabel="P1–P6 axes" />
        <StatTile
          label="Hotspot modules"
          value={data.n_hotspot_modules.toLocaleString("en-US")}
          sublabel="per-cell-type"
        />
        <StatTile
          label="Novel modules"
          value={data.n_novel.toLocaleString("en-US")}
          sublabel="vs known programs"
        />
        <StatTile
          label="cNMF change-points"
          value={data.changepoints.length}
          sublabel="unbiased detection"
        />
      </div>

      {/* Framing note */}
      <div className="rounded-lg border border-dashed border-border bg-muted/20 p-4 text-sm text-muted-foreground">
        The k=6 NMF programs are presented as{" "}
        <span className="font-medium text-foreground">
          continuous interpretive axes (P1–P6)
        </span>
        , not discrete patient subtypes (2026-07-04 audit): the factor subspace is
        cross-seed stable, but no internal metric selects k=6 as a hard partition.
        Loading bars show each program&rsquo;s curated defining genes (confounder-aware).
      </div>

      {/* Program × stage activity heatmap */}
      <section>
        <h2 className="mb-1 text-lg font-semibold">Program activity across fibrosis stage</h2>
        <p className="mb-3 max-w-2xl text-sm text-muted-foreground">
          Mean NMF program score per fibrosis stage (F0–F4). Darker = higher mean
          activity across samples at that stage.
        </p>
        <div className="rounded-lg border border-border bg-card p-4">
          <Heatmap
            data={heatmapCells}
            rows={programs.map((p) => p.code)}
            cols={STAGES}
            colorScale="sequential"
            showValues
            valueFormat={(v) => v.toFixed(3)}
            ariaLabel="Program by fibrosis stage activity heatmap"
          />
          <div className="mt-2 flex flex-wrap gap-x-4 gap-y-1">
            {programs.map((p) => (
              <span
                key={p.code}
                className="inline-flex items-center gap-1.5 text-xs text-muted-foreground"
              >
                <span
                  aria-hidden
                  className="size-2.5 rounded-full"
                  style={{ backgroundColor: programColor(p.code) }}
                />
                <span className="font-medium text-foreground">{p.code}</span> {p.label}
              </span>
            ))}
          </div>
        </div>
      </section>

      {/* Program cards */}
      <section>
        <h2 className="mb-3 text-lg font-semibold">Programs &amp; defining genes</h2>
        <div className="grid gap-5 md:grid-cols-2">
          {programs.map((p) => (
            <ProgramCard
              key={p.code}
              program={p}
              bars={topGenes[p.code] ?? []}
              scoreDomain={scoreDomain}
            />
          ))}
        </div>
      </section>

      {/* Hotspot modules */}
      <section>
        <h2 className="mb-1 text-lg font-semibold">
          Hotspot autocorrelation modules
        </h2>
        <p className="mb-3 max-w-2xl text-sm text-muted-foreground">
          Per-cell-type Hotspot modules (DeTomaso &amp; Yosef 2021). {data.n_novel} of{" "}
          {data.n_hotspot_modules} are novel vs cNMF / bulk-NMF / Hallmark / SCENIC+
          references.
        </p>
        <DataTable
          data={data.hotspot_modules}
          columns={hotspotCols}
          pageSize={15}
          initialSort={{ key: "stability_score", dir: "desc" }}
          rowKey={(r) => `${r.cell_type}:${r.module}`}
        />
      </section>

      {/* Change-points */}
      <section>
        <h2 className="mb-1 text-lg font-semibold">cNMF program change-points</h2>
        <p className="mb-3 max-w-2xl text-sm text-muted-foreground">
          Unbiased change-point detection on cNMF global k=16 programs (not the k=6
          NMF axes above): best model per program and the stage at which activity
          shifts, where a segmented fit wins.
        </p>
        <DataTable
          data={data.changepoints}
          columns={changepointCols}
          pageSize={8}
          initialSort={{ key: "delta_aic", dir: "asc" }}
          rowKey={(r) => r.program}
        />
      </section>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Sex-biased tab
// ---------------------------------------------------------------------------

function SexTab({ data }: { data: SexSummary }) {
  // Shared gene order: strongest sex-stratum effect first.
  const ordered = useMemo(() => {
    const displayLabel = (g: SexGene) =>
      g.symbol && g.symbol !== "nan" ? g.symbol : `Unnamed·${g.chr}`;
    return [...data.genes]
      .map((g) => ({ g, label: displayLabel(g) }))
      .sort(
        (a, b) =>
          Math.max(Math.abs(b.g.logfc_male), Math.abs(b.g.logfc_female)) -
          Math.max(Math.abs(a.g.logfc_male), Math.abs(a.g.logfc_female))
      );
  }, [data.genes]);

  const femaleBars: DivergingBarRow[] = ordered.map(({ g, label }) => ({
    label,
    value: g.logfc_female,
    sig: g.interaction_padj < 0.05,
  }));
  const maleBars: DivergingBarRow[] = ordered.map(({ g, label }) => ({
    label,
    value: g.logfc_male,
    sig: g.interaction_padj < 0.05,
  }));

  const tableCols: DataTableColumn<SexGene>[] = [
    {
      key: "symbol",
      header: "Gene",
      sortable: true,
      render: (r) =>
        r.symbol && r.symbol !== "nan" ? (
          <GeneLink symbol={r.symbol} />
        ) : (
          <span className="font-mono text-sm italic text-muted-foreground">
            Unnamed
          </span>
        ),
    },
    {
      key: "sex_class",
      header: "Sex bias",
      sortable: true,
      render: (r) => (
        <ProgramChip
          label={prettySexClass(r.sex_class)}
          color={SEX_COLORS[r.sex_class] ?? CONTROL}
        />
      ),
    },
    {
      key: "biotype",
      header: "Biotype",
      render: (r) => (
        <span className="text-xs text-muted-foreground">{r.biotype}</span>
      ),
    },
    {
      key: "chr",
      header: "Chr",
      render: (r) => <span className="font-mono text-xs">{r.chr}</span>,
    },
    {
      key: "logfc_female",
      header: "logFC (F)",
      numeric: true,
      sortable: true,
      render: (r) => fmtSignedLFC(r.logfc_female),
    },
    {
      key: "logfc_male",
      header: "logFC (M)",
      numeric: true,
      sortable: true,
      render: (r) => fmtSignedLFC(r.logfc_male),
    },
    {
      key: "interaction_padj",
      header: "Interaction padj",
      numeric: true,
      sortable: true,
      render: (r) => fmtInteractionPadj(r.interaction_padj),
    },
  ];

  return (
    <div className="space-y-8">
      <div className="grid grid-cols-3 gap-4">
        <StatTile label="Sex-dimorphic genes" value={data.n_dimorphic} />
        <StatTile
          label="Male-biased"
          value={data.n_male_biased}
          sublabel="larger effect in males"
        />
        <StatTile
          label="Female-biased"
          value={data.n_female_biased}
          sublabel="larger effect in females"
        />
      </div>

      <div className="rounded-lg border border-dashed border-border bg-muted/20 p-4 text-sm text-muted-foreground">
        <span className="font-medium text-foreground">Sex-biased framing only.</span>{" "}
        {data.note}
      </div>

      <section>
        <div className="mb-3 flex flex-wrap items-center justify-between gap-2">
          <h2 className="text-lg font-semibold">
            Female vs male stratum effect per gene
          </h2>
          <Legend
            items={[
              { color: DIVERGING.down, label: "down (− logFC)" },
              { color: CONTROL, label: "zero (control)" },
              { color: DIVERGING.up, label: "up (+ logFC)" },
            ]}
          />
        </div>
        <div className="grid gap-6 md:grid-cols-2">
          <div className="rounded-lg border border-border bg-card p-4">
            <DivergingBar
              data={femaleBars}
              title="Female stratum"
              xLabel="log₂FC (female)"
              clamp={2}
              showValues
              ariaLabel="Female-stratum log fold-change per dimorphic gene"
            />
          </div>
          <div className="rounded-lg border border-border bg-card p-4">
            <DivergingBar
              data={maleBars}
              title="Male stratum"
              xLabel="log₂FC (male)"
              clamp={2}
              showValues
              ariaLabel="Male-stratum log fold-change per dimorphic gene"
            />
          </div>
        </div>
      </section>

      <section>
        <h2 className="mb-3 text-lg font-semibold">Dimorphic gene detail</h2>
        <DataTable
          data={data.genes}
          columns={tableCols}
          pageSize={0}
          initialSort={{ key: "interaction_padj", dir: "asc" }}
          rowKey={(r, i) => `${r.symbol}:${r.chr}:${i}`}
        />
      </section>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Spatial observability tab
// ---------------------------------------------------------------------------

function SpatialProgramsTab({
  coverage,
  effects,
  unavailable,
}: {
  coverage: SpatialProgramCoverageRow[];
  effects: SpatialProgramEffectRow[];
  unavailable: boolean;
}) {
  const programOptions = useMemo(
    () =>
      Array.from(
        new Map(
          coverage.map((row) => [
            row.program_uid,
            { uid: row.program_uid, label: row.program_label },
          ])
        ).values()
      ).sort((a, b) => a.label.localeCompare(b.label)),
    [coverage]
  );
  const [programUid, setProgramUid] = useState("");
  const selectedUid = programOptions.some((row) => row.uid === programUid)
    ? programUid
    : programOptions[0]?.uid ?? "";
  const selectedCoverage = coverage.filter((row) => row.program_uid === selectedUid);
  const selectedEffects = effects.filter((row) => row.program_uid === selectedUid);

  if (unavailable || coverage.length === 0) {
    return (
      <EmptyState
        title="Spatial Resource tables are not staged"
        description="The portal consumer is ready, but categorical spatial tables remain in the isolated candidate until synchronized Resource promotion."
      />
    );
  }

  const coverageCols: DataTableColumn<SpatialProgramCoverageRow>[] = [
    { key: "dataset_id", header: "Dataset", sortable: true },
    {
      key: "coverage_status",
      header: "Coverage",
      sortable: true,
      render: (r) => <ProgramChip label={r.coverage_status.replace(/_/g, " ")} />,
    },
    { key: "n_genes_measured", header: "Genes", numeric: true, sortable: true },
    {
      key: "retained_l1_weight",
      header: "Retained L1",
      numeric: true,
      sortable: true,
      render: (r) => r.retained_l1_weight.toFixed(3),
    },
    {
      key: "biological_unit_resolution",
      header: "Unit resolution",
      sortable: true,
      render: (r) => r.biological_unit_resolution.replace(/_/g, " "),
    },
    {
      key: "dataset_gate",
      header: "Dataset gate",
      sortable: true,
      render: (r) => r.dataset_gate.replace(/_/g, " "),
    },
    {
      key: "testability_reason",
      header: "Reason",
      render: (r) => r.testability_reason.replace(/_/g, " "),
    },
  ];
  const effectCols: DataTableColumn<SpatialProgramEffectRow>[] = [
    { key: "dataset_id", header: "Dataset", sortable: true },
    {
      key: "evidence_state",
      header: "Evidence state",
      sortable: true,
      render: (r) => <ProgramChip label={r.evidence_state.replace(/_/g, " ")} />,
    },
    {
      key: "estimate",
      header: "Estimate",
      numeric: true,
      sortable: true,
      render: (r) => (r.estimate == null ? "—" : r.estimate.toFixed(4)),
    },
    { key: "effect_unit", header: "Effect unit" },
    {
      key: "qvalue",
      header: "q",
      numeric: true,
      sortable: true,
      render: (r) => (r.qvalue == null ? "—" : r.qvalue.toPrecision(3)),
    },
    {
      key: "n_biological",
      header: "Biological n",
      numeric: true,
      render: (r) => r.n_biological ?? "—",
    },
    {
      key: "n_technical",
      header: "Technical n",
      numeric: true,
      render: (r) => r.n_technical ?? "—",
    },
  ];

  return (
    <div className="space-y-8">
      <div className="rounded-lg border border-dashed border-border bg-muted/20 p-4 text-sm text-muted-foreground">
        All 117 frozen programs receive outcome-free coverage. Confirmatory inference is
        restricted to the two predeclared robust-display programs. Indeterminate and
        untestable states are retained; Moran&apos;s I is spatial organization, not disease direction.
      </div>
      <label className="block max-w-xl text-sm font-medium">
        Frozen program
        <select
          className="mt-1 w-full rounded-md border border-border bg-background px-3 py-2 text-sm"
          value={selectedUid}
          onChange={(event) => setProgramUid(event.target.value)}
        >
          {programOptions.map((row) => (
            <option key={row.uid} value={row.uid}>
              {row.label} · {row.uid}
            </option>
          ))}
        </select>
      </label>
      <section>
        <h2 className="mb-1 text-lg font-semibold">Assay observability</h2>
        <p className="mb-3 text-sm text-muted-foreground">
          Release {selectedCoverage[0]?.release_id}. Biological-unit resolution and source gates
          determine whether coverage can support inference.
        </p>
        <DataTable
          data={selectedCoverage}
          columns={coverageCols}
          pageSize={0}
          initialSort={{ key: "dataset_id", dir: "asc" }}
          rowKey={(r) => `${r.dataset_id}:${r.assay_id}:${r.program_uid}`}
        />
      </section>
      <section>
        <h2 className="mb-1 text-lg font-semibold">Confirmatory spatial effects</h2>
        {selectedEffects.length > 0 ? (
          <DataTable
            data={selectedEffects}
            columns={effectCols}
            pageSize={0}
            initialSort={{ key: "dataset_id", dir: "asc" }}
            rowKey={(r) => `${r.dataset_id}:${r.assay_id}:${r.program_uid}`}
          />
        ) : (
          <p className="text-sm text-muted-foreground">
            Coverage only. This program was not part of the sealed two-program confirmatory family.
          </p>
        )}
      </section>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Root client
// ---------------------------------------------------------------------------

const TABS = [
  { id: "programs", label: "Programs" },
  { id: "spatial", label: "Spatial observability" },
  { id: "sex", label: "Sex-biased" },
];

export function ProgramsClient() {
  const [tab, setTab] = useState<string>("programs");
  const [programs, setPrograms] = useState<ProgramsSummary | null>(null);
  const [sex, setSex] = useState<SexSummary | null>(null);
  const [topGenes, setTopGenes] = useState<Record<string, MembershipRow[]>>({});
  const [spatialCoverage, setSpatialCoverage] = useState<SpatialProgramCoverageRow[]>([]);
  const [spatialEffects, setSpatialEffects] = useState<SpatialProgramEffectRow[]>([]);
  const [spatialUnavailable, setSpatialUnavailable] = useState(false);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(false);

  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        const [prog, sexData] = await Promise.all([
          fetch(dataUrl("programs_summary.json")).then((r) => {
            if (!r.ok) throw new Error(`programs_summary.json ${r.status}`);
            return r.json() as Promise<ProgramsSummary>;
          }),
          fetch(dataUrl("sex_summary.json")).then((r) => {
            if (!r.ok) throw new Error(`sex_summary.json ${r.status}`);
            return r.json() as Promise<SexSummary>;
          }),
        ]);
        if (cancelled) return;
        setPrograms(prog);
        setSex(sexData);

        // Per-program loading bars: restrict to each program's curated top_genes
        // (confounder-aware) and read their per-program logFC from the parquet.
        const curatedByProgram = new Map<string, Set<string>>();
        for (const p of prog.nmf_programs) {
          curatedByProgram.set(p.code, new Set(p.top_genes.slice(0, 14)));
        }
        const union = Array.from(
          new Set(prog.nmf_programs.flatMap((p) => p.top_genes.slice(0, 14)))
        );
        const grouped: Record<string, MembershipRow[]> = {};
        if (union.length > 0) {
          const inList = union.map((g) => sqlString(g)).join(",");
          const rows = await queryParquet<MembershipRow>(
            "program_gene_membership.parquet",
            `SELECT program_code, symbol, logfc, padj, direction
             FROM program_gene_membership
             WHERE symbol IN (${inList})`
          );
          if (cancelled) return;
          for (const r of rows) {
            const allowed = curatedByProgram.get(r.program_code);
            if (!allowed || !allowed.has(r.symbol)) continue;
            (grouped[r.program_code] ??= []).push(r);
          }
          for (const code of Object.keys(grouped)) {
            grouped[code] = grouped[code]
              .sort((a, b) => Math.abs(b.logfc) - Math.abs(a.logfc))
              .slice(0, TOP_GENES_PER_PROGRAM);
          }
        }
        if (cancelled) return;
        setTopGenes(grouped);
      } catch {
        if (!cancelled) setError(true);
      } finally {
        if (!cancelled) setLoading(false);
      }
    })();
    return () => {
      cancelled = true;
    };
  }, []);

  useEffect(() => {
    let cancelled = false;
    Promise.all([
      queryParquet<SpatialProgramCoverageRow>(
        "spatial_program_coverage.parquet",
        `SELECT release_id, dataset_id, assay_id, program_uid, program_label,
                n_program_genes, n_genes_measured, retained_l1_weight,
                coverage_status, testability_reason, dataset_gate, biological_unit,
                biological_unit_resolution, source_dependence
           FROM spatial_program_coverage`
      ),
      queryParquet<SpatialProgramEffectRow>(
        "spatial_program_effects.parquet",
        `SELECT release_id, dataset_id, assay_id, program_uid, estimand, effect_unit,
                estimate, qvalue, n_biological, n_technical, evidence_state,
                testability_reason
           FROM spatial_program_effects`
      ),
    ])
      .then(([coverageRows, effectRows]) => {
        if (cancelled) return;
        setSpatialCoverage(coverageRows);
        setSpatialEffects(effectRows);
      })
      .catch(() => {
        if (!cancelled) setSpatialUnavailable(true);
      });
    return () => {
      cancelled = true;
    };
  }, []);

  return (
    <div className="mt-6">
      <TabBar
        tabs={TABS}
        active={tab}
        onChange={setTab}
        aria-label="Programs sections"
        className="mb-6"
      />

      {loading ? (
        <div className="space-y-4">
          <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
            {Array.from({ length: 4 }).map((_, i) => (
              <SkeletonBlock key={i} className="h-20" />
            ))}
          </div>
          <SkeletonBlock className="h-56" />
          <SkeletonBlock className="h-72" />
        </div>
      ) : error || (!programs && !sex) ? (
        <EmptyState
          title="Could not load program data"
          description="programs_summary.json / sex_summary.json failed to load. Check the data directory."
        />
      ) : tab === "programs" ? (
        programs ? (
          <ProgramsTab data={programs} topGenes={topGenes} />
        ) : (
          <EmptyState title="No program data" />
        )
      ) : tab === "spatial" ? (
        <SpatialProgramsTab
          coverage={spatialCoverage}
          effects={spatialEffects}
          unavailable={spatialUnavailable}
        />
      ) : sex ? (
        <SexTab data={sex} />
      ) : (
        <EmptyState title="No sex-biased data" />
      )}
    </div>
  );
}
