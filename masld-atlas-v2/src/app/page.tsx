"use client";

/**
 * Landing page — the "Living Atlas" composition.
 *
 * Hero (gene constellation) → animated KPI band → full-bleed gene ticker →
 * bento entry grid of real-data mini-viz previews → at-a-glance data story →
 * featured genes. Each section is a self-contained client component; the heavy
 * previews lazy-mount inside their bento tiles, so first paint stays cheap.
 */

import { Hero } from "@/components/hero/hero";
import { HeroStats } from "@/components/landing/hero-stats";
import { GeneTicker } from "@/components/landing/gene-ticker";
import { BentoGrid } from "@/components/landing/bento-grid";
import { AtlasStory } from "@/components/landing/atlas-story";
import { FeaturedGenes } from "@/components/landing/featured-genes";

export default function HomePage() {
  return (
    <div className="pb-16">
      <Hero />

      {/* Animated KPI band, tucked up into the hero's lower vignette. */}
      <div className="mx-auto max-w-6xl px-6">
        <HeroStats className="relative z-10 -mt-8 sm:-mt-12" />
      </div>

      {/* Full-bleed gene ticker — top differentially expressed genes. */}
      <div className="mt-10 border-y border-border/60 bg-card/40">
        <GeneTicker className="py-3" />
      </div>

      <div className="mx-auto max-w-6xl px-6 pt-12">
        {/* Bento entry grid. */}
        <section className="mb-12">
          <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-muted-foreground">
            Explore the atlas
          </h2>
          <BentoGrid />
        </section>

        <AtlasStory />

        <FeaturedGenes />
      </div>
    </div>
  );
}
