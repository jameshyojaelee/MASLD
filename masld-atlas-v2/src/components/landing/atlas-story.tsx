"use client";

/**
 * AtlasStory — a compact scroll-reveal narrative strip below the bento grid.
 *
 * A staggered fade-rise of "beats", each a big mono number + one honest
 * caption. Every number is imported from atlas-constants (never a hardcoded
 * digit). An optional leading "N samples across M cohorts" beat is added only
 * when `atlas_summary.json` resolves (fail-soft; no hardcoded sample/cohort
 * count). Ends with a CTA into the atlas.
 */

import { useEffect, useState } from "react";
import { HashLink as Link } from "@/components/hash-link";
import { Reveal, RevealItem } from "@/components/motion/reveal";
import { dataUrl } from "@/lib/data-base";
import { cn } from "@/lib/utils";
import {
  ATLAS_GENES,
  DEG_COUNT,
  DEG_GATE_LABEL,
  CONVERGENCE_TIER1,
  DRUGS_APPROVED,
  DRUGS_CLINICAL,
  DRUGS_PRECLINICAL,
  fmt,
} from "@/lib/atlas-constants";
import type { AtlasSummary } from "@/lib/types";

interface Beat {
  value: string;
  caption: string;
}

export function AtlasStory() {
  const [summary, setSummary] = useState<AtlasSummary | null>(null);

  useEffect(() => {
    let alive = true;
    fetch(dataUrl("atlas_summary.json"))
      .then((r) => r.json())
      .then((d: AtlasSummary) => {
        if (alive) setSummary(d);
      })
      .catch(() => {});
    return () => {
      alive = false;
    };
  }, []);

  const beats: Beat[] = [];

  if (summary?.total_samples && summary?.total_cohorts) {
    beats.push({
      value: fmt(summary.total_samples),
      caption: `samples across ${summary.total_cohorts} control-bearing cohorts`,
    });
  }

  beats.push(
    { value: fmt(ATLAS_GENES), caption: "genes profiled across the atlas" },
    { value: fmt(DEG_COUNT), caption: `differential genes (${DEG_GATE_LABEL})` },
    {
      value: fmt(CONVERGENCE_TIER1),
      caption: "Tier-1 convergent targets with multi-omics agreement",
    },
    {
      value: String(DRUGS_APPROVED),
      caption: `approved drugs recovered (${fmt(DRUGS_CLINICAL)} clinical · ${fmt(
        DRUGS_PRECLINICAL
      )} preclinical)`,
    }
  );

  const colClass =
    beats.length >= 5
      ? "sm:grid-cols-3 lg:grid-cols-5"
      : "sm:grid-cols-2 lg:grid-cols-4";

  return (
    <section className="mb-12">
      <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-muted-foreground">
        At a glance
      </h2>
      <Reveal className={cn("grid gap-3", colClass)}>
        {beats.map((beat) => (
          <RevealItem
            key={beat.caption}
            className="rounded-lg border border-border bg-card p-4"
          >
            <div className="font-mono text-3xl font-semibold tracking-tight tabular-nums">
              {beat.value}
            </div>
            <p className="mt-1 text-xs text-muted-foreground leading-relaxed">
              {beat.caption}
            </p>
          </RevealItem>
        ))}
      </Reveal>
      <Reveal stagger={false} className="mt-4 flex justify-end">
        <Link
          href="/atlas"
          className="text-sm font-medium text-primary hover:underline"
        >
          Explore the atlas &rarr;
        </Link>
      </Reveal>
    </section>
  );
}
