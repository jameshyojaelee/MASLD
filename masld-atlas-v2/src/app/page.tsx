"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Separator } from "@/components/ui/separator";
import { EvidenceFingerprint } from "@/components/evidence-fingerprint";
import { useAppStore } from "@/lib/store";
import type { FeaturedGene } from "@/lib/types";

const QUICK_LINKS = [
  {
    title: "Atlas",
    description: "10-cohort Integrated mega-analysis with 5,484 DEGs",
    href: "/atlas",
    color: "text-chart-1",
  },
  {
    title: "Progression",
    description: "F0-F4 disease trajectory and the F2 switch",
    href: "/progression",
    color: "text-chart-2",
  },
  {
    title: "Causal Architecture",
    description: "COLOC + TWAS across 24 GWAS studies",
    href: "/causal",
    color: "text-chart-3",
  },
  {
    title: "Drug Pipeline",
    description: "133 reversal compounds, 15 validated targets",
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
  const { setCommandOpen } = useAppStore();
  const [featured, setFeatured] = useState<FeaturedGene[]>([]);

  useEffect(() => {
    fetch("/data/featured_genes.json")
      .then((r) => r.json())
      .then(setFeatured)
      .catch(() => {});
  }, []);

  return (
    <div className="mx-auto max-w-5xl px-6 py-10">
      {/* Hero */}
      <section className="mb-12 text-center">
        <h1 className="text-4xl font-bold tracking-tight sm:text-5xl">
          MASLD Atlas
        </h1>
        <p className="mt-4 text-lg text-muted-foreground">
          Multi-modal atlas for metabolic dysfunction-associated steatotic
          liver disease.
        </p>
        <div className="mt-6 flex items-center justify-center gap-3">
          <Button size="lg" onClick={() => setCommandOpen(true)}>
            <svg className="size-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M21 21l-5.197-5.197m0 0A7.5 7.5 0 105.196 5.196a7.5 7.5 0 0010.607 10.607z" />
            </svg>
            Search genes
          </Button>
          <Button variant="outline" size="lg" render={<Link href="/atlas" />}>
            Explore atlas
          </Button>
        </div>
      </section>

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
                href={`/gene/${gene.symbol}/`}
                className="group rounded-lg border border-border bg-card p-4 transition-colors hover:border-primary/40"
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
              className="rounded-lg border border-border bg-card p-4 transition-colors hover:border-primary/40"
            >
              <p className={`text-sm font-semibold ${link.color}`}>{link.title}</p>
              <p className="mt-1 text-xs text-muted-foreground">{link.description}</p>
            </Link>
          ))}
        </div>
      </section>

    </div>
  );
}
