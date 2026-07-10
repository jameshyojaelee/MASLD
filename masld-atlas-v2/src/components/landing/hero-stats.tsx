"use client";

/**
 * <HeroStats> — the landing KPI row.
 *
 * A responsive grid (2 → 3 → 6 columns) of hero <StatTile>s, each animated in on
 * scroll via <Reveal>/<RevealItem> (stagger). Every displayed number comes from
 * `atlas-constants` — never hardcode a digit here so a results refresh is a
 * one-file edit. Numeric tiles CountUp; the tiny approved-drugs integer renders
 * statically. Icons are muted UI chrome; no text is colored.
 */

import { Database, Dna, Layers, Microscope, Pill, Target } from "lucide-react";
import { Reveal, RevealItem } from "@/components/motion/reveal";
import { StatTile } from "@/components/stat-tile";
import {
  ATLAS_GENES,
  COLOC_SUSIE,
  COLOC_UNION,
  CONVERGENCE_TIER1,
  DEG_COUNT,
  DEG_DOWN,
  DEG_UP,
  DRUGS_APPROVED,
  DRUGS_CLINICAL,
  DRUGS_PRECLINICAL,
  GWAS_COUNT,
  SC_CELLS,
  fmt,
} from "@/lib/atlas-constants";

/** Compact millions formatter, e.g. 1232318 -> "1.23M". */
const fmtMillions = (n: number): string => `${(n / 1e6).toFixed(2)}M`;

export function HeroStats({ className }: { className?: string }) {
  return (
    <Reveal className={className}>
      <div className="grid grid-cols-2 gap-3 sm:grid-cols-3 lg:grid-cols-6">
        <RevealItem>
          <StatTile
            variant="hero"
            label="DEGs"
            numericValue={DEG_COUNT}
            sublabel={`↑${fmt(DEG_UP)} / ↓${fmt(DEG_DOWN)}`}
            icon={<Dna />}
          />
        </RevealItem>
        <RevealItem>
          <StatTile
            variant="hero"
            label="COLOC effectors"
            numericValue={COLOC_SUSIE}
            sublabel={`${fmt(COLOC_UNION)} union · ${GWAS_COUNT} GWAS`}
            icon={<Target />}
          />
        </RevealItem>
        <RevealItem>
          <StatTile
            variant="hero"
            label="Convergent targets"
            numericValue={CONVERGENCE_TIER1}
            sublabel="Tier-1 multi-omics"
            icon={<Layers />}
          />
        </RevealItem>
        <RevealItem>
          <StatTile
            variant="hero"
            label="Approved drugs"
            value={DRUGS_APPROVED}
            sublabel={`${fmt(DRUGS_CLINICAL)} clinical · ${fmt(
              DRUGS_PRECLINICAL
            )} preclinical`}
            icon={<Pill />}
          />
        </RevealItem>
        <RevealItem>
          <StatTile
            variant="hero"
            label="Genes profiled"
            numericValue={ATLAS_GENES}
            icon={<Database />}
          />
        </RevealItem>
        <RevealItem>
          <StatTile
            variant="hero"
            label="Cells"
            numericValue={SC_CELLS}
            format={fmtMillions}
            icon={<Microscope />}
          />
        </RevealItem>
      </div>
    </Reveal>
  );
}

export default HeroStats;
