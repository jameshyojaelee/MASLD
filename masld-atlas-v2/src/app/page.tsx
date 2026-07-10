"use client";

import { useEffect, useState } from "react";
import { HashLink as Link } from "@/components/hash-link";
import { Badge } from "@/components/ui/badge";
import { Separator } from "@/components/ui/separator";
import { EvidenceFingerprint } from "@/components/evidence-fingerprint";
import { Hero } from "@/components/hero/hero";
import { dataUrl } from "@/lib/data-base";
import {
  DEG_COUNT,
  GWAS_COUNT,
  DRUGS_APPROVED,
  DRUGS_CLINICAL,
  DRUGS_PRECLINICAL,
  fmt,
} from "@/lib/atlas-constants";
import type { FeaturedGene } from "@/lib/types";

const QUICK_LINKS = [
  {
    title: "Atlas",
    description: `Pooled 5-cohort analysis · ${fmt(DEG_COUNT)} DEGs`,
    href: "/atlas",
    color: "text-chart-1",
  },
  {
    title: "Progression",
    description: "F0-F4 disease trajectory and stage transitions",
    href: "/progression",
    color: "text-chart-2",
  },
  {
    title: "Causal Architecture",
    description: `COLOC + TWAS across ${GWAS_COUNT} GWAS (5 ancestries)`,
    href: "/genetics",
    color: "text-chart-3",
  },
  {
    title: "Drug Pipeline",
    description: `${DRUGS_APPROVED} approved · ${fmt(DRUGS_CLINICAL)} clinical · ${fmt(DRUGS_PRECLINICAL)} preclinical`,
    href: "/drugs",
    color: "text-chart-4",
  },
  {
    title: "Gene Explorer",
    description: "Search and compare genes across 7 evidence layers",
    href: "/explore",
    color: "text-chart-5",
  },
  {
    title: "Downloads",
    description: "Full atlas, gene lists, and supplementary tables",
    href: "/downloads",
    color: "text-muted-foreground",
  },
];

export default function HomePage() {
  const [featured, setFeatured] = useState<FeaturedGene[]>([]);

  useEffect(() => {
    fetch(dataUrl("featured_genes.json"))
      .then((r) => r.json())
      .then(setFeatured)
      .catch(() => {});
  }, []);

  return (
    <div className="pb-16">
      <Hero />

      <div className="mx-auto max-w-6xl px-6">
        <Separator className="mb-12" />

        {/* Featured genes */}
        {featured.length > 0 && (
          <section className="mb-12">
            <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-muted-foreground">
              Featured genes
            </h2>
            <div className="grid gap-4 sm:grid-cols-2 lg:grid-cols-3">
              {featured.map((gene) => (
                <Link
                  key={gene.symbol}
                  href={`#/gene?symbol=${gene.symbol}`}
                  className="hover-lift group rounded-lg border border-border bg-card p-4 hover:border-primary/40"
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
                  <span>logFC: <span className="font-mono">{gene.bulk_logfc.toFixed(3)}</span></span>
                  {gene.coloc_pp4 != null && (
                    <span>PP.H4: <span className="font-mono">{gene.coloc_pp4.toFixed(3)}</span></span>
                  )}
                  {gene.drug && (
                    <span className="truncate">{gene.drug}</span>
                  )}
                </div>
              </Link>
            ))}
          </div>
        </section>
      )}

      <Separator className="mb-12" />

      {/* Quick links */}
      <section className="mb-12">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-muted-foreground">
          Quick links
        </h2>
        <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-3">
          {QUICK_LINKS.map((link) => (
            <Link
              key={link.href}
              href={link.href}
              className="hover-lift rounded-lg border border-border bg-card p-4 hover:border-primary/40"
            >
              <p className={`text-sm font-semibold ${link.color}`}>{link.title}</p>
              <p className="mt-1 text-xs text-muted-foreground">{link.description}</p>
            </Link>
          ))}
        </div>
      </section>
      </div>
    </div>
  );
}
