"use client";

import { useEffect, useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";
import { PageContainer } from "@/components/page-container";
import { PageHeader } from "@/components/page-header";
import { Volcano, type VolcanoPoint } from "@/components/charts";
import { SkeletonBlock } from "@/components/states";
import { getGeneIndex } from "@/lib/search-index";
import { ATLAS_GENES, DEG_COUNT, DEG_GATE_LABEL, fmt } from "@/lib/atlas-constants";
import type { GeneIndexEntry } from "@/lib/types";

// ---------------------------------------------------------------------------
// Hardcoded study metadata
// ---------------------------------------------------------------------------

interface Cohort {
  dataset: string;
  samples: number;
  condition: string;
  platform: string;
  reference: string;
}

const COHORTS: Cohort[] = [
  {
    dataset: "GSE135251",
    samples: 216,
    condition: "MASLD staging (NAS)",
    platform: "RNA-seq",
    reference: "Govaere 2020",
  },
  {
    dataset: "GSE130970",
    samples: 78,
    condition: "NASH vs Healthy",
    platform: "RNA-seq",
    reference: "Hoang 2019",
  },
  {
    dataset: "GSE126848",
    samples: 112,
    condition: "NAFL/NASH/Fibrosis",
    platform: "RNA-seq",
    reference: "Suppli 2019",
  },
  {
    dataset: "GSE167523",
    samples: 118,
    condition: "NAFLD progression",
    platform: "RNA-seq",
    reference: "Pantano 2021",
  },
  {
    dataset: "GSE162694",
    samples: 76,
    condition: "Steatosis/NASH",
    platform: "RNA-seq",
    reference: "Fang 2021",
  },
  {
    dataset: "GSE174478",
    samples: 64,
    condition: "NAFL/NASH",
    platform: "RNA-seq",
    reference: "Vvedenskaya 2021",
  },
  {
    dataset: "GSE193066",
    samples: 216,
    condition: "MASLD staging",
    platform: "RNA-seq",
    reference: "Govaere 2022",
  },
  {
    dataset: "GSE213621",
    samples: 256,
    condition: "Fibrosis staging",
    platform: "RNA-seq",
    reference: "Vali 2023",
  },
  {
    dataset: "GSE240729",
    samples: 108,
    condition: "Fibrosis",
    platform: "RNA-seq",
    reference: "Liu 2024",
  },
];

// PRJNA512027 (Gerhard 2018) collected but dropped from presentation for an
// L0/S0 library-prep batch confound; retained in the pipeline config for
// Script 05d provenance only, not shown here.
const TOTAL_SAMPLES = 1244;

// ---------------------------------------------------------------------------
// Pipeline steps
// ---------------------------------------------------------------------------

interface PipelineStep {
  label: string;
  description: string;
}

const PIPELINE_STEPS: PipelineStep[] = [
  {
    label: "Per-study DE",
    description: "limma-voom per cohort with study-specific contrasts",
  },
  {
    label: "QC Filtering",
    description: "PCA + library size filters (pass_technical)",
  },
  {
    label: "Count Integration",
    description: "Harmonized count matrix across 9 cohorts",
  },
  {
    label: "Pooled (Cohort-Adjusted) Analysis",
    description: "limma-voom quality-weighted, dataset as a fixed effect",
  },
  {
    label: "DEG Thresholding",
    description: DEG_GATE_LABEL,
  },
  {
    label: "Consensus DEGs",
    description: `${fmt(DEG_COUNT)} DEGs, leave-one-out cross-validated`,
  },
];

// ---------------------------------------------------------------------------
// Key result cards
// ---------------------------------------------------------------------------

interface ResultCard {
  value: string;
  label: string;
  sublabel: string;
}

const RESULT_CARDS: ResultCard[] = [
  {
    value: fmt(ATLAS_GENES),
    label: "Genes tested",
    sublabel: "GENCODE v49 / GRCh38",
  },
  {
    value: fmt(DEG_COUNT),
    label: "DEGs",
    sublabel: "fdr < 0.05 at lfc = 0.25",
  },
  {
    value: "81.6%",
    label: "LOO-CV recovery",
    sublabel: "rho = 0.952 · 5-fold leave-one-cohort-out",
  },
  {
    value: "7,842",
    label: "Robust genes",
    sublabel: "Significant in all 5 of 5 LOO folds",
  },
];

// ---------------------------------------------------------------------------
// Page
// ---------------------------------------------------------------------------

export default function AtlasPage() {
  const router = useRouter();
  const [genes, setGenes] = useState<GeneIndexEntry[]>([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    getGeneIndex().then((data) => {
      setGenes(data);
      setLoading(false);
    });
  }, []);

  // Build the volcano dataset. Significance is the canonical DEG flag (the
  // effect-size-aware interval-null gate), passed explicitly so the shared
  // Volcano colors up = red / down = blue / non-sig = control gray.
  const volcanoData = useMemo<VolcanoPoint[]>(() => {
    const out: VolcanoPoint[] = [];
    for (const g of genes) {
      if (g.bulk_logfc == null || g.bulk_padj == null) continue;
      out.push({
        symbol: g.symbol,
        logFC: g.bulk_logfc,
        padj: g.bulk_padj,
        sig: g.is_deg,
      });
    }
    return out;
  }, [genes]);

  return (
    <PageContainer>
      <PageHeader
        title="Atlas Construction"
        description={
          <>
            5-cohort pooled (cohort-adjusted) analysis across 846 samples —{" "}
            {fmt(DEG_COUNT)} DEGs by the {DEG_GATE_LABEL}, validated by
            leave-one-out cross-validation.
          </>
        }
      />

      {/* ------------------------------------------------------------------ */}
      {/* Section 1: Study Overview                                            */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Study Overview
        </h2>
        <div className="overflow-x-auto rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50">
                {["Dataset", "Samples", "Condition", "Platform", "Reference"].map(
                  (h) => (
                    <th
                      key={h}
                      className="px-4 py-2.5 text-left text-xs font-medium text-muted-foreground"
                    >
                      {h}
                    </th>
                  )
                )}
              </tr>
            </thead>
            <tbody>
              {COHORTS.map((c) => (
                <tr
                  key={c.dataset}
                  className="border-b border-border/50 transition-colors hover:bg-muted/30"
                >
                  <td className="px-4 py-2.5 font-mono text-xs font-semibold text-primary">
                    {c.dataset}
                  </td>
                  <td className="px-4 py-2.5 text-right font-mono text-xs">
                    {c.samples}
                  </td>
                  <td className="px-4 py-2.5 text-xs">{c.condition}</td>
                  <td className="px-4 py-2.5">
                    <Badge variant="outline" className="text-[10px]">
                      {c.platform}
                    </Badge>
                  </td>
                  <td className="px-4 py-2.5 text-xs text-muted-foreground">
                    {c.reference}
                  </td>
                </tr>
              ))}
              {/* Total row */}
              <tr className="bg-muted/30">
                <td className="px-4 py-2.5 text-xs font-semibold">
                  Total (QC-passing)
                </td>
                <td className="px-4 py-2.5 text-right font-mono text-xs font-semibold">
                  {TOTAL_SAMPLES.toLocaleString()}
                </td>
                <td colSpan={3} className="px-4 py-2.5 text-xs text-muted-foreground">
                  9 independent cohorts
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 2: Analysis Pipeline                                         */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Analysis Pipeline
        </h2>
        {/* Horizontal step flow */}
        <div className="flex flex-wrap items-start gap-0">
          {PIPELINE_STEPS.map((step, idx) => (
            <div key={step.label} className="flex items-start">
              {/* Step card */}
              <div className="flex w-[148px] flex-col gap-1.5 rounded-lg border border-border bg-card px-3 py-3 shadow-sm">
                <span className="text-[10px] font-semibold uppercase tracking-wider text-muted-foreground">
                  Step {idx + 1}
                </span>
                <span className="text-sm font-semibold leading-tight">
                  {step.label}
                </span>
                <span className="text-[11px] leading-tight text-muted-foreground">
                  {step.description}
                </span>
              </div>
              {/* Arrow connector */}
              {idx < PIPELINE_STEPS.length - 1 && (
                <div className="flex h-[72px] items-center px-1 text-muted-foreground">
                  <span className="text-lg">&#8594;</span>
                </div>
              )}
            </div>
          ))}
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 3: Key Results                                               */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Key Results
        </h2>
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
          {RESULT_CARDS.map((card) => (
            <div
              key={card.label}
              className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm hover-lift"
            >
              <p className="font-numeric text-2xl font-bold text-primary">
                {card.value}
              </p>
              <p className="mt-0.5 text-sm font-semibold">{card.label}</p>
              <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                {card.sublabel}
              </p>
            </div>
          ))}
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 4: Volcano Plot                                              */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <div className="mb-3 flex items-center justify-between">
          <h2 className="text-xl font-semibold tracking-tight">Volcano Plot</h2>
          {!loading && (
            <span className="text-xs text-muted-foreground">
              {volcanoData.length.toLocaleString()} genes plotted
            </span>
          )}
        </div>
        <div className="rounded-lg border border-border bg-card p-4 shadow-sm">
          {loading ? (
            <div className="space-y-3">
              <SkeletonBlock className="h-[420px] w-full rounded-lg" />
              <div className="flex gap-4">
                <SkeletonBlock className="h-3 w-28" />
                <SkeletonBlock className="h-3 w-32" />
                <SkeletonBlock className="h-3 w-24" />
              </div>
            </div>
          ) : (
            <Volcano
              data={volcanoData}
              height={460}
              ariaLabel="Volcano plot of disease vs healthy differential expression"
              onPointClick={(sym) =>
                router.push(`/gene?symbol=${encodeURIComponent(sym)}`)
              }
              caption={
                <>
                  Direction by hue (red = up, blue = down) with significance by
                  opacity; non-significant genes fall back to control gray.
                  Dashed guides mark |logFC|&nbsp;=&nbsp;0.25 and
                  FDR&nbsp;=&nbsp;0.05. Click any colored point to open the gene
                  profile. Points at the top edge represent
                  padj&nbsp;&lt;&nbsp;1e-300 (clamped for display).
                </>
              }
            />
          )}
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Footer nav                                                           */}
      {/* ------------------------------------------------------------------ */}
      <div className="flex flex-wrap gap-3 border-t border-border pt-6">
        <Link
          href="/explore"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Browse DEGs in Explorer &rarr;
        </Link>
        <Link
          href="/causal"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Causal Architecture &rarr;
        </Link>
      </div>
    </PageContainer>
  );
}
