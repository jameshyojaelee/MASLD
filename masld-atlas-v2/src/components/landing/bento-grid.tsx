"use client";

/**
 * The landing entry grid: an asymmetric bento of the atlas's primary sections.
 *
 * Replaces the flat quick-links row. The Single-cell UMAP tile is a 2x2 hero;
 * the rest are 1x1 cells with tiny real-data previews (ancestry / funnel /
 * volcano / stage), a ⌘K gene-search prompt, and two nav tiles with faint
 * decorative watermarks. Layout
 * is a plain CSS grid (1 col mobile → 2 → 4) with explicit spans for rhythm.
 * Section stats come from `atlas-constants`; previews are lazy-mounted inside
 * each tile so they never block first paint.
 */

import { Search, ArrowRight, Download, ArrowLeftRight } from "lucide-react";
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
      {/* Hero: Single-cell UMAP (2×2) — the flagship atlas figure */}
      <BentoTile
        title="Single-cell"
        description="Integrated liver cell atlas across 16 cell types"
        stat={`${scMillions}M cells`}
        href="/single-cell"
        className="[background-image:var(--gradient-surface)] sm:col-span-2 sm:row-span-2 lg:col-span-2 lg:row-span-2"
      >
        <MiniUmap />
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

      {/* Atlas — bulk differential expression */}
      <BentoTile
        title="Atlas"
        description="Pooled 5-cohort disease-vs-control differential expression"
        stat={`${fmt(DEG_COUNT)} DEGs`}
        href="/atlas"
      >
        <MiniVolcano />
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

      {/* Downloads — no data preview; a faint watermark icon gives it weight */}
      <BentoTile
        title="Downloads"
        description="Full atlas, gene lists, supplementary tables"
        stat={<ArrowRight className="size-3.5" aria-hidden />}
        href="/downloads"
        lazy={false}
      >
        <div className="relative h-full w-full">
          <Download
            aria-hidden
            strokeWidth={1.25}
            className="pointer-events-none absolute -bottom-3 -right-2 size-24 text-muted-foreground/15"
          />
        </div>
      </BentoTile>

      {/* Translation — cross-species arrows as a faint decorative watermark */}
      <BentoTile
        title="Translation"
        description="Cross-species and clinical translation"
        stat={<ArrowRight className="size-3.5" aria-hidden />}
        href="/translation"
        lazy={false}
      >
        <div className="relative h-full w-full">
          <ArrowLeftRight
            aria-hidden
            strokeWidth={1.25}
            className="pointer-events-none absolute -bottom-3 -right-2 size-24 text-muted-foreground/15"
          />
        </div>
      </BentoTile>
    </div>
  );
}
