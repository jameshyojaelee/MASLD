import Link from "next/link";
import { Badge } from "@/components/ui/badge";

// ---------------------------------------------------------------------------
// Deconvolution stat cards
// ---------------------------------------------------------------------------

interface StatCard {
  value: string;
  label: string;
  sublabel: string;
}

const DECONV_CARDS: StatCard[] = [
  {
    value: "1,918",
    label: "Significant DEGs tested",
    sublabel: "MuSiC deconvolution across 1,444 samples",
  },
  {
    value: "81.5%",
    label: "Hepatocyte-intrinsic",
    sublabel: "1,563 / 1,918 significant DEGs",
  },
  {
    value: "18.5%",
    label: "Immune / stromal-attributed",
    sublabel: "355 / 1,918 significant DEGs",
  },
  {
    value: "0.931",
    label: "BayesPrism concordance",
    sublabel: "Jaccard index between methods",
  },
];

// ---------------------------------------------------------------------------
// SCENIC+ regulon table
// ---------------------------------------------------------------------------

interface RegulonRow {
  tf: string;
  regulon_size: number;
  activity_diff: number;
  direction: "Up in MASLD" | "Down in MASLD";
}

const REGULON_TABLE: RegulonRow[] = [
  { tf: "HNF4A", regulon_size: 487, activity_diff: -0.15, direction: "Down in MASLD" },
  { tf: "RORA", regulon_size: 312, activity_diff: -0.12, direction: "Down in MASLD" },
  { tf: "THRB", regulon_size: 205, activity_diff: -0.08, direction: "Down in MASLD" },
  { tf: "CEBPA", regulon_size: 156, activity_diff: -0.11, direction: "Down in MASLD" },
  { tf: "NR1H4", regulon_size: 98, activity_diff: -0.09, direction: "Down in MASLD" },
  { tf: "STAT3", regulon_size: 234, activity_diff: 0.18, direction: "Up in MASLD" },
  { tf: "NFE2L2", regulon_size: 189, activity_diff: 0.14, direction: "Up in MASLD" },
  { tf: "JUN", regulon_size: 167, activity_diff: 0.16, direction: "Up in MASLD" },
];

// ---------------------------------------------------------------------------
// Sex stratification stat cards
// ---------------------------------------------------------------------------

const SEX_CARDS: StatCard[] = [
  {
    value: "1,978",
    label: "Female-biased DEGs",
    sublabel: "Integrated interaction model (group × sex), padj < 0.05",
  },
  {
    value: "675",
    label: "Male-biased DEGs",
    sublabel: "Integrated interaction model (group × sex), padj < 0.05",
  },
  {
    value: "427",
    label: "Divergent DEGs",
    sublabel: "Opposite direction between sexes",
  },
  {
    value: "2.9:1",
    label: "Female : Male ratio",
    sublabel: "3,080 total sex-dimorphic (interaction padj < 0.05)",
  },
];

// ---------------------------------------------------------------------------
// SVG: Donut chart for deconvolution
// ---------------------------------------------------------------------------

function DeconvDonut() {
  const size = 200;
  const cx = size / 2;
  const cy = size / 2;
  const outerR = 85;
  const innerR = 55;
  const hepatocyteFrac = 0.815;

  // Compute arc for hepatocyte segment (81.5%)
  const hepatoAngle = hepatocyteFrac * 2 * Math.PI;
  const startAngle = -Math.PI / 2; // start at top

  // Hepatocyte arc end
  const hepatoEndAngle = startAngle + hepatoAngle;
  const hepatoX1 = cx + outerR * Math.cos(startAngle);
  const hepatoY1 = cy + outerR * Math.sin(startAngle);
  const hepatoX2 = cx + outerR * Math.cos(hepatoEndAngle);
  const hepatoY2 = cy + outerR * Math.sin(hepatoEndAngle);
  const hepatoIX1 = cx + innerR * Math.cos(hepatoEndAngle);
  const hepatoIY1 = cy + innerR * Math.sin(hepatoEndAngle);
  const hepatoIX2 = cx + innerR * Math.cos(startAngle);
  const hepatoIY2 = cy + innerR * Math.sin(startAngle);

  // Other arc (18.5%)
  const otherEndAngle = hepatoEndAngle + (1 - hepatocyteFrac) * 2 * Math.PI;
  const otherX2 = cx + outerR * Math.cos(otherEndAngle);
  const otherY2 = cy + outerR * Math.sin(otherEndAngle);
  const otherIX1 = cx + innerR * Math.cos(otherEndAngle);
  const otherIY1 = cy + innerR * Math.sin(otherEndAngle);

  const hepatoPath = [
    `M ${hepatoX1} ${hepatoY1}`,
    `A ${outerR} ${outerR} 0 1 1 ${hepatoX2} ${hepatoY2}`,
    `L ${hepatoIX1} ${hepatoIY1}`,
    `A ${innerR} ${innerR} 0 1 0 ${hepatoIX2} ${hepatoIY2}`,
    "Z",
  ].join(" ");

  const otherPath = [
    `M ${hepatoX2} ${hepatoY2}`,
    `A ${outerR} ${outerR} 0 0 1 ${otherX2} ${otherY2}`,
    `L ${otherIX1} ${otherIY1}`,
    `A ${innerR} ${innerR} 0 0 0 ${hepatoIX1} ${hepatoIY1}`,
    "Z",
  ].join(" ");

  return (
    <svg
      width={size}
      height={size}
      viewBox={`0 0 ${size} ${size}`}
      className="shrink-0"
      role="img"
      aria-label="Donut chart: 81.5% hepatocyte-intrinsic, 18.5% immune/stromal"
    >
      <path d={hepatoPath} fill="hsl(210 100% 56%)" opacity={0.85} />
      <path d={otherPath} fill="hsl(var(--muted-foreground))" opacity={0.4} />
      <text
        x={cx}
        y={cy - 6}
        textAnchor="middle"
        className="fill-foreground text-[22px] font-bold"
      >
        81.5%
      </text>
      <text
        x={cx}
        y={cy + 12}
        textAnchor="middle"
        className="fill-muted-foreground text-[10px]"
      >
        Hepatocyte
      </text>
    </svg>
  );
}

// ---------------------------------------------------------------------------
// SVG: Bar chart for sex stratification
// ---------------------------------------------------------------------------

function SexBarChart() {
  const svgW = 500;
  const svgH = 180;
  const padLeft = 100;
  const padRight = 60;
  const padTop = 20;
  const padBottom = 40;
  const innerW = svgW - padLeft - padRight;
  const barH = 28;
  const gap = 16;

  const data = [
    { label: "Female-biased", count: 1978, color: "#f472b6" },
    { label: "Divergent", count: 427, color: "#a1a1aa" },
    { label: "Male-biased", count: 675, color: "#60a5fa" },
  ];

  const maxCount = Math.max(...data.map((d) => d.count));

  return (
    <svg
      width={svgW}
      height={svgH}
      viewBox={`0 0 ${svgW} ${svgH}`}
      className="w-full max-w-xl"
      role="img"
      aria-label="Bar chart: 1,978 female-biased, 427 divergent, 675 male-biased DEGs (interaction model)"
    >
      {data.map((d, i) => {
        const y = padTop + i * (barH + gap);
        const w = Math.max(4, (d.count / maxCount) * innerW);
        return (
          <g key={d.label}>
            <text
              x={padLeft - 8}
              y={y + barH / 2 + 4}
              textAnchor="end"
              className="fill-muted-foreground text-[11px]"
            >
              {d.label}
            </text>
            <rect
              x={padLeft}
              y={y}
              width={w}
              height={barH}
              rx={4}
              fill={d.color}
              opacity={0.8}
            />
            <text
              x={padLeft + w + 8}
              y={y + barH / 2 + 4}
              className="fill-foreground text-[12px] font-semibold"
            >
              {d.count.toLocaleString()}
            </text>
          </g>
        );
      })}
      {/* X axis */}
      <line
        x1={padLeft}
        x2={svgW - padRight}
        y1={svgH - padBottom}
        y2={svgH - padBottom}
        stroke="hsl(var(--border))"
        strokeWidth={1}
      />
      <text
        x={padLeft + innerW / 2}
        y={svgH - 10}
        textAnchor="middle"
        className="fill-muted-foreground text-[11px]"
      >
        Number of DEGs
      </text>
    </svg>
  );
}

// ---------------------------------------------------------------------------
// Page
// ---------------------------------------------------------------------------

export default function ResolutionPage() {
  return (
    <div className="mx-auto max-w-5xl px-6 py-10">
      {/* ------------------------------------------------------------------ */}
      {/* Header                                                               */}
      {/* ------------------------------------------------------------------ */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Cell-Type &amp; Sex Resolution
        </h1>
        <p className="mt-2 text-muted-foreground">
          Deconvolution attribution, SCENIC+ regulon analysis, and sex-stratified
          Integrated mega-analysis reveal cell-type-specific and sex-dimorphic disease
          signatures.
        </p>
      </div>

      {/* ------------------------------------------------------------------ */}
      {/* Section 1: Deconvolution Attribution                                 */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Deconvolution Attribution
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          MuSiC deconvolution across 1,444 samples reveals that 81.5%
          (1,563/1,918) of significant DEGs are hepatocyte-intrinsic. BayesPrism
          concordance Jaccard&nbsp;=&nbsp;0.931.
        </p>

        {/* Stat cards */}
        <div className="mb-6 grid grid-cols-2 gap-4 sm:grid-cols-4">
          {DECONV_CARDS.map((card) => (
            <div
              key={card.label}
              className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm"
            >
              <p className="font-mono text-2xl font-bold text-primary">
                {card.value}
              </p>
              <p className="mt-0.5 text-sm font-semibold">{card.label}</p>
              <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                {card.sublabel}
              </p>
            </div>
          ))}
        </div>

        {/* Donut chart */}
        <div className="flex items-center gap-6 rounded-lg border border-border bg-muted/30 p-6">
          <DeconvDonut />
          <div className="space-y-2 text-sm">
            <div className="flex items-center gap-2">
              <span
                className="inline-block h-3 w-3 rounded-sm"
                style={{ background: "hsl(210 100% 56%)", opacity: 0.85 }}
              />
              <span>
                <span className="font-semibold">Hepatocyte-intrinsic</span>{" "}
                <span className="text-muted-foreground">
                  &mdash; 1,563 DEGs (81.5%)
                </span>
              </span>
            </div>
            <div className="flex items-center gap-2">
              <span
                className="inline-block h-3 w-3 rounded-sm"
                style={{
                  background: "hsl(var(--muted-foreground))",
                  opacity: 0.4,
                }}
              />
              <span>
                <span className="font-semibold">Immune / stromal</span>{" "}
                <span className="text-muted-foreground">
                  &mdash; 355 DEGs (18.5%)
                </span>
              </span>
            </div>
            <p className="mt-2 text-xs text-muted-foreground">
              Two independent methods (MuSiC + BayesPrism) yield Jaccard =
              0.931 concordance on hepatocyte-intrinsic classification.
            </p>
          </div>
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 2: SCENIC+ Regulon Analysis                                  */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          SCENIC+ Regulon Analysis
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          139 hepatocyte regulons identified via SCENIC+ from scATAC-seq. 25
          disease-associated regulons show differential activity between MASLD and
          healthy.
        </p>

        <div className="overflow-x-auto rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50">
                {["TF", "Regulon Size", "Activity Diff", "Direction"].map(
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
              {REGULON_TABLE.map((row) => {
                const isUp = row.direction === "Up in MASLD";
                return (
                  <tr
                    key={row.tf}
                    className="border-b border-border/50 transition-colors hover:bg-muted/30"
                  >
                    <td className="px-4 py-2.5">
                      <Link
                        href={`/gene/${encodeURIComponent(row.tf)}`}
                        className="font-mono text-xs font-semibold text-primary hover:underline"
                      >
                        {row.tf}
                      </Link>
                    </td>
                    <td className="px-4 py-2.5 text-right font-mono text-xs">
                      {row.regulon_size} targets
                    </td>
                    <td className="px-4 py-2.5 text-right font-mono text-xs">
                      <span
                        className={
                          isUp
                            ? "text-red-500 dark:text-red-400"
                            : "text-blue-500 dark:text-blue-400"
                        }
                      >
                        {row.activity_diff >= 0 ? "+" : ""}
                        {row.activity_diff.toFixed(2)}
                      </span>
                    </td>
                    <td className="px-4 py-2.5">
                      <Badge
                        variant={isUp ? "destructive" : "secondary"}
                        className="text-[10px]"
                      >
                        {row.direction}
                      </Badge>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
        <p className="mt-2 text-xs text-muted-foreground">
          Top 8 disease-associated regulons by absolute activity difference.
          Click any TF name to view its gene profile with full evidence.
        </p>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 3: Sex Stratification                                        */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Sex Stratification
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          Sex-stratified Integrated interaction model (group_binary × inferred_sex)
          across matched-dataset samples reveals 3,080 sex-dimorphic DEGs (interaction padj &lt; 0.05).
        </p>

        {/* Stat cards */}
        <div className="mb-6 grid grid-cols-2 gap-4 sm:grid-cols-4">
          {SEX_CARDS.map((card) => (
            <div
              key={card.label}
              className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm"
            >
              <p className="font-mono text-2xl font-bold text-primary">
                {card.value}
              </p>
              <p className="mt-0.5 text-sm font-semibold">{card.label}</p>
              <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                {card.sublabel}
              </p>
            </div>
          ))}
        </div>

        {/* Bar chart */}
        <div className="mb-4 rounded-lg border border-border bg-muted/30 p-6">
          <SexBarChart />
          <div className="mt-3 flex gap-4 text-xs text-muted-foreground">
            <span className="flex items-center gap-1.5">
              <span className="inline-block h-2.5 w-2.5 rounded-full bg-pink-400" />
              Female-biased (1,978)
            </span>
            <span className="flex items-center gap-1.5">
              <span className="inline-block h-2.5 w-2.5 rounded-full bg-zinc-400" />
              Divergent (427)
            </span>
            <span className="flex items-center gap-1.5">
              <span className="inline-block h-2.5 w-2.5 rounded-full bg-blue-400" />
              Male-biased (675)
            </span>
          </div>
        </div>

        {/* Key finding callout */}
        <div className="rounded-lg border border-amber-300 bg-amber-50 px-5 py-4 dark:border-amber-700 dark:bg-amber-950/30">
          <p className="text-sm font-semibold text-amber-900 dark:text-amber-200">
            Key Finding
          </p>
          <p className="mt-1 text-sm text-amber-800 dark:text-amber-300">
            Female-biased DEGs are <strong>depleted</strong> for COLOC
            (OR&nbsp;=&nbsp;0.64, padj&nbsp;=&nbsp;0.001), indicating female
            transcriptomic effects are non-genetic. Male-biased DEGs show no
            such depletion (OR&nbsp;=&nbsp;0.94, n.s.).
          </p>
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Footer nav                                                           */}
      {/* ------------------------------------------------------------------ */}
      <div className="flex flex-wrap gap-3 border-t border-border pt-6">
        <Link
          href="/atlas"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          &larr; Atlas Construction
        </Link>
        <Link
          href="/causal"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Causal Architecture &rarr;
        </Link>
      </div>
    </div>
  );
}
