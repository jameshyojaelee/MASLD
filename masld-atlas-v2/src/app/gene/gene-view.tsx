"use client";

import { useRef, useState } from "react";
import { useReducedMotion } from "framer-motion";
import { EvidenceBar } from "@/components/charts";
import { StatTile } from "@/components/stat-tile";
import { TabBar, type TabItem } from "@/components/tab-bar";
import {
  ProgramChip,
  TierBadge,
} from "@/components/chips";
import { EmptyState, SkeletonBlock } from "@/components/states";
import { DEG_GATE_LABEL } from "@/lib/atlas-constants";
import { divergingColor, CONTROL, type ModalityKey } from "@/lib/palette";
import {
  useParquetRows,
  evidenceStrengths,
  num,
  fmtLogFC,
  fmtP,
  fmtPP4,
  directionArrow,
  DASH,
  type AtlasCoreRow,
} from "./use-gene-data";
import {
  ExpressionSection,
  GeneticsSection,
  SingleCellSection,
  SpatialSection,
  ProgramsSection,
  TherapeuticsSection,
} from "./gene-sections";
import { sqlString } from "@/lib/duck";

const TABS: TabItem[] = [
  { id: "expression", label: "Expression" },
  { id: "genetics", label: "Genetics" },
  { id: "singlecell", label: "Single-cell" },
  { id: "spatial", label: "Spatial" },
  { id: "programs", label: "Programs" },
  { id: "therapeutics", label: "Therapeutics" },
];

// Evidence-fingerprint bar → the detail tab that best covers that modality.
// (Not every modality has a dedicated tab; these route to the closest section.)
const MODALITY_TAB: Record<ModalityKey, string> = {
  s1_human: "expression",
  s2_genetic: "genetics",
  s3_essential: "therapeutics",
  s4_epigenomic: "genetics",
  s5_spatial: "spatial",
  s6_singlecell: "singlecell",
  s7_mouse: "expression",
  s8_proteomics: "expression",
};

export function GeneView({ symbol }: { symbol: string }) {
  const [tab, setTab] = useState("expression");
  const tabsRef = useRef<HTMLDivElement>(null);
  const reduceMotion = useReducedMotion();

  const goToModality = (key: ModalityKey) => {
    setTab(MODALITY_TAB[key] ?? "expression");
    requestAnimationFrame(() =>
      tabsRef.current?.scrollIntoView({
        behavior: reduceMotion ? "auto" : "smooth",
        block: "start",
      })
    );
  };

  const header = useParquetRows<AtlasCoreRow>(
    "atlas_core.parquet",
    `SELECT * FROM atlas_core WHERE upper(human_symbol) = upper(${sqlString(
      symbol
    )}) LIMIT 1`
  );

  if (header.error) {
    return (
      <div className="mx-auto max-w-6xl px-6 py-10">
        <EmptyState
          title="Could not load gene data"
          description={header.error}
        />
      </div>
    );
  }

  if (header.loading) {
    return (
      <div className="mx-auto max-w-6xl space-y-4 px-6 py-8">
        <SkeletonBlock className="h-8 w-48" />
        <SkeletonBlock className="h-24" />
        <SkeletonBlock className="h-64" />
      </div>
    );
  }

  const row = header.data?.[0];
  if (!row) {
    return (
      <div className="mx-auto max-w-6xl px-6 py-10">
        <EmptyState
          title={`No gene “${symbol}” in the atlas`}
          description="Check the symbol, or search from the command palette (⌘K)."
        />
      </div>
    );
  }

  const logfc = num(row.bulk_logFC);
  const bestPP4 = Math.max(
    num(row.coloc_best_susie_pp4) ?? 0,
    num(row.coloc_abf_best_pp4) ?? 0
  );
  const degColor = row.is_deg && logfc != null ? divergingColor(logfc) : CONTROL;

  return (
    <div className="mx-auto max-w-6xl px-6 py-8">
      {/* ---------------------------------------------------------------- */}
      {/* Identity + status                                                */}
      {/* ---------------------------------------------------------------- */}
      <div className="mb-6">
        <div className="flex flex-wrap items-baseline gap-x-3 gap-y-1">
          <h1 className="font-mono text-3xl font-bold italic tracking-tight">
            {row.human_symbol}
          </h1>
          {row.gene_biotype && (
            <span className="text-sm text-muted-foreground">
              {row.gene_biotype.replace(/_/g, " ")}
            </span>
          )}
          {row.ensembl_id && (
            <span className="font-numeric text-xs text-muted-foreground">
              {row.ensembl_id}
            </span>
          )}
        </div>

        <div className="mt-3 flex flex-wrap items-center gap-2">
          <ProgramChip
            color={degColor}
            label={
              row.is_deg
                ? `${directionArrow(logfc)} Differentially expressed`
                : "Not a DEG"
            }
          />
          {row.convergence_tier && <TierBadge tier={row.convergence_tier} />}
          {row.concordance_state && (
            <ProgramChip label={row.concordance_state.replace(/_/g, " ")} />
          )}
          {row.coloc_cross_ancestry_replicated && (
            <ProgramChip label="Cross-ancestry replicated" />
          )}
          {row.drug_dev_status && row.drug_dev_status !== "none" && (
            <ProgramChip label={`Drug: ${row.drug_dev_status}`} />
          )}
        </div>

        <p className="mt-2 text-xs text-muted-foreground">
          DEG calls use the {DEG_GATE_LABEL}.
        </p>
      </div>

      {/* ---------------------------------------------------------------- */}
      {/* Headline stats + evidence fingerprint                            */}
      {/* ---------------------------------------------------------------- */}
      <div className="mb-6 grid gap-4 lg:grid-cols-[2fr_1fr]">
        <div className="grid grid-cols-2 gap-3 sm:grid-cols-3">
          <StatTile
            label="Bulk log2 FC"
            value={fmtLogFC(logfc)}
            sublabel="disease vs control"
            delta={
              logfc != null
                ? {
                    value: directionArrow(logfc),
                    direction: logfc > 0 ? "up" : logfc < 0 ? "down" : "neutral",
                  }
                : undefined
            }
          />
          <StatTile
            label="FDR"
            value={fmtP(row.treat_fdr)}
            sublabel="interval-null gate"
          />
          <StatTile
            label="Coloc PP.H4"
            value={bestPP4 > 0 ? fmtPP4(bestPP4) : DASH}
            sublabel={row.coloc_best_susie_gwas ?? "no colocalization"}
          />
          <StatTile
            label="Convergence rank"
            value={num(row.convergence_rank) != null ? String(num(row.convergence_rank)) : DASH}
            sublabel={
              row.convergence_score != null
                ? `score ${fmtLogFC(row.convergence_score)}`
                : undefined
            }
          />
          <StatTile
            label="Essentiality"
            value={fmtLogFC(row.essentiality_chronos)}
            sublabel="DepMap CHRONOS"
          />
          <StatTile
            label="Active layers"
            value={num(row.layers_active) != null ? String(num(row.layers_active)) : DASH}
            sublabel="evidence modalities"
          />
        </div>

        <div className="rounded-lg border border-border bg-card p-4 shadow-sm">
          <h3 className="mb-3 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
            Evidence fingerprint
          </h3>
          <EvidenceBar
            strengths={evidenceStrengths(row)}
            showLabels
            onModalityClick={goToModality}
          />
        </div>
      </div>

      {/* ---------------------------------------------------------------- */}
      {/* Detail tabs (each queries its own long-format parquet lazily)    */}
      {/* ---------------------------------------------------------------- */}
      <div ref={tabsRef} className="scroll-mt-20">
        <TabBar
          tabs={TABS}
          active={tab}
          onChange={setTab}
          aria-label="Gene evidence sections"
          className="mb-4"
        />
      </div>

      {tab === "expression" && <ExpressionSection symbol={row.human_symbol} />}
      {tab === "genetics" && <GeneticsSection symbol={row.human_symbol} />}
      {tab === "singlecell" && <SingleCellSection symbol={row.human_symbol} />}
      {tab === "spatial" && <SpatialSection symbol={row.human_symbol} />}
      {tab === "programs" && <ProgramsSection symbol={row.human_symbol} />}
      {tab === "therapeutics" && <TherapeuticsSection symbol={row.human_symbol} />}
    </div>
  );
}
