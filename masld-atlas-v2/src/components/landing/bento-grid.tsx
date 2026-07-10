"use client";

/**
 * The landing entry grid: an asymmetric bento of the atlas's primary sections.
 *
 * Replaces the flat quick-links row. The Atlas tile is a 2x2 hero; the rest are
 * 1x1 cells with tiny real-data previews (volcano / ancestry / funnel / UMAP /
 * stage), a ⌘K gene-search prompt, and two text-only navigation tiles. Layout
 * is a plain CSS grid (1 col mobile → 2 → 4) with explicit spans for rhythm.
 * Section stats come from `atlas-constants`; previews are lazy-mounted inside
 * each tile so they never block first paint.
 */

import { Search, ArrowRight } from "lucide-react";
import { BentoTile } from "./bento-tile";
import { MiniVolcano } from "@/components/previews/mini-volcano";
import { AncestryBars } from "@/components/previews/ancestry-bars";
import { DrugFunnel } from "@/components/previews/drug-funnel";
import { MiniUmap } from "@/components/previews/mini-umap";
import { StageStrip } from "@/components/previews/stage-strip";
import { useAppStore } from "@/lib/store";
import {
  DEG_COUNT,
  COLOC_SUSIE,
  GWAS_COUNT,
  DRUGS_APPROVED,
  SC_CELLS,
  fmt,
} from "@/lib/atlas-constants";

export function BentoGrid() {
  const setCommandOpen = useAppStore((s) => s.setCommandOpen);
  const scMillions = (SC_CELLS / 1e6).toFixed(2);

  return (
    <div className="grid auto-rows-[172px] grid-cols-1 gap-3 sm:grid-cols-2 lg:grid-cols-4">
      {/* Hero: Atlas (2×2) */}
      <BentoTile
        title="Atlas"
        description="Pooled 5-cohort disease-vs-control differential expression"
        stat={`${fmt(DEG_COUNT)} DEGs`}
        href="/atlas"
        className="[background-image:var(--gradient-surface)] sm:col-span-2 sm:row-span-2 lg:col-span-2 lg:row-span-2"
      >
        <MiniVolcano />
      </BentoTile>

      {/* Genetics */}
      <BentoTile
        title="Genetics"
        description="COLOC + TWAS causal architecture across ancestries"
        stat={`${COLOC_SUSIE} effectors · ${GWAS_COUNT} GWAS`}
        href="/genetics"
      >
        <AncestryBars />
      </BentoTile>

      {/* Drugs */}
      <BentoTile
        title="Drugs"
        description="Target-to-drug development pipeline"
        stat={`${DRUGS_APPROVED} approved`}
        href="/drugs"
      >
        <DrugFunnel />
      </BentoTile>

      {/* Single-cell */}
      <BentoTile
        title="Single-cell"
        description="Integrated liver cell atlas"
        stat={`${scMillions}M cells`}
        href="/single-cell"
      >
        <MiniUmap />
      </BentoTile>

      {/* Progression */}
      <BentoTile
        title="Progression"
        description="Fibrosis-stage disease trajectory"
        stat="F0–F4"
        href="/progression"
      >
        <StageStrip />
      </BentoTile>

      {/* Gene Explorer — ⌘K search prompt (2×1) */}
      <BentoTile
        title="Gene Explorer"
        description="Search any gene across 7 evidence layers"
        onClick={() => setCommandOpen(true)}
        lazy={false}
        className="sm:col-span-2 lg:col-span-2"
      >
        <div className="flex h-full items-center">
          <div className="flex w-full items-center gap-2 rounded-md border border-border bg-background/60 px-3 py-2 text-sm text-muted-foreground transition-colors group-hover:border-primary/40">
            <Search className="size-4 shrink-0" aria-hidden />
            <span className="truncate">Search genes, e.g. HKDC1, THBS2…</span>
            <kbd className="ml-auto shrink-0 rounded border border-border bg-muted px-1.5 py-0.5 font-numeric text-[10px] text-muted-foreground">
              ⌘K
            </kbd>
          </div>
        </div>
      </BentoTile>

      {/* Downloads (text only) */}
      <BentoTile
        title="Downloads"
        description="Full atlas, gene lists, supplementary tables"
        stat={<ArrowRight className="size-3.5" aria-hidden />}
        href="/downloads"
      />

      {/* Translation (text only) */}
      <BentoTile
        title="Translation"
        description="Cross-species and clinical translation"
        stat={<ArrowRight className="size-3.5" aria-hidden />}
        href="/translation"
      />
    </div>
  );
}
