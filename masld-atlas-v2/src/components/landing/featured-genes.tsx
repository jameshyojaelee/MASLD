"use client";

/**
 * FeaturedGenes — landing-page strip of highlighted atlas genes.
 *
 * Fetches `featured_genes.json` and renders a scroll-revealed grid of cards
 * (mono symbol + category badge + evidence fingerprint + tagline + logFC, plus
 * a lead drug when the gene carries one). Cards deep-link to the gene view.
 * Fail-soft: renders nothing while loading, on fetch error, or when empty.
 *
 * Extracted from the old inline block in page.tsx, which read `gene.coloc_pp4`
 * and `gene.drug` — fields absent from the JSON (always undefined). Those dead
 * reads are gone; a real target now comes from the `drugs[]` array.
 */

import { useEffect, useState } from "react";
import { HashLink as Link } from "@/components/hash-link";
import { Badge } from "@/components/ui/badge";
import { EvidenceFingerprint } from "@/components/evidence-fingerprint";
import { Reveal, RevealItem } from "@/components/motion/reveal";
import { dataUrl } from "@/lib/data-base";
import type { FeaturedGene } from "@/lib/types";

export function FeaturedGenes() {
  const [featured, setFeatured] = useState<FeaturedGene[]>([]);

  useEffect(() => {
    let alive = true;
    fetch(dataUrl("featured_genes.json"))
      .then((r) => r.json())
      .then((d: FeaturedGene[]) => {
        if (alive) setFeatured(Array.isArray(d) ? d : []);
      })
      .catch(() => {});
    return () => {
      alive = false;
    };
  }, []);

  if (featured.length === 0) return null;

  return (
    <section className="mb-12">
      <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-muted-foreground">
        Featured genes
      </h2>
      <Reveal className="grid gap-4 sm:grid-cols-2 lg:grid-cols-3">
        {featured.map((gene) => {
          const lead = gene.drugs?.[0];
          return (
            <RevealItem key={gene.symbol}>
              <Link
                href={`#/gene?symbol=${gene.symbol}`}
                className="hover-lift group block h-full rounded-lg border border-border bg-card p-4 hover:border-primary/40"
              >
                <div className="flex items-start justify-between">
                  <div>
                    <span className="font-mono text-base font-semibold group-hover:text-primary">
                      {gene.symbol}
                    </span>
                    <Badge variant="secondary" className="ml-2">
                      {gene.category}
                    </Badge>
                  </div>
                  <EvidenceFingerprint evidence={gene.evidence} size={36} showTooltip={false} />
                </div>
                <p className="mt-2 text-xs text-muted-foreground leading-relaxed">
                  {gene.tagline}
                </p>
                <div className="mt-3 flex items-center gap-3 text-xs text-muted-foreground">
                  <span>
                    logFC: <span className="font-mono">{gene.bulk_logfc.toFixed(3)}</span>
                  </span>
                  {lead && (
                    <span className="truncate">
                      {lead.drug}
                      <span className="text-muted-foreground/70"> · {lead.stage}</span>
                    </span>
                  )}
                </div>
              </Link>
            </RevealItem>
          );
        })}
      </Reveal>
    </section>
  );
}
