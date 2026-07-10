import { HashLink as Link } from "@/components/hash-link";
import { Badge } from "@/components/ui/badge";
import { PageContainer } from "@/components/page-container";
import { PageHeader } from "@/components/page-header";
import { categoricalColor, CONTROL } from "@/lib/palette";

// LEGACY: this route is absorbed by /single-cell — retire after migration. Kept
// live for now; changes here are correctness-only (no expanded surface).

// Qualitative two-class palette (not directional): hepatocyte vs other, and the
// male / female sex split. All routed through the palette authority.
const HEPATO_COLOR = categoricalColor(0); // blue
const MALE_COLOR = categoricalColor(0); // blue
const FEMALE_COLOR = categoricalColor(3); // reddish purple

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
    value: "1,497",
    label: "Classifiable DEGs",
    sublabel: "MuSiC deconvolution, pooled (cohort-adjusted) analysis",
  },
  {
    value: "65.1%",
    label: "Hepatocyte-intrinsic",
    sublabel: "975 / 1,497 classifiable DEGs",
  },
  {
    value: "34.9%",
    label: "Immune / stromal-attributed",
    sublabel: "522 / 1,497 classifiable DEGs",
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
    value: "7",
    label: "Male-biased DEGs",
    sublabel: "LVQW-fixed (C2) interaction model, padj < 0.05",
  },
  {
    value: "1",
    label: "Female-biased DEGs",
    sublabel: "LVQW-fixed (C2) interaction model, padj < 0.05",
  },
  {
    value: "8",
    label: "Total sex-dimorphic",
    sublabel: "Down from 290 after fixing the dataset term (C2)",
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
  const hepatocyteFrac = 0.651;

  // Compute arc for hepatocyte segment (65.1%)
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

  // Other arc (34.9%)
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
      aria-label="Donut chart: 65.1% hepatocyte-intrinsic, 34.9% immune/stromal"
    >
      <path d={hepatoPath} fill={HEPATO_COLOR} opacity={0.85} />
      <path d={otherPath} fill={CONTROL} opacity={0.55} />
      <text
        x={cx}
        y={cy - 6}
        textAnchor="middle"
        className="fill-foreground text-[22px] font-bold"
      >
        65.1%
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
    { label: "Male-biased", count: 7, color: MALE_COLOR },
    { label: "Female-biased", count: 1, color: FEMALE_COLOR },
  ];

  const maxCount = Math.max(...data.map((d) => d.count));

  return (
    <svg
      width={svgW}
      height={svgH}
      viewBox={`0 0 ${svgW} ${svgH}`}
      className="w-full max-w-xl"
      role="img"
      aria-label="Bar chart: 7 male-biased, 1 female-biased DEGs (LVQW-fixed interaction model)"
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
        stroke="var(--color-border)"
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
    <PageContainer>
      {/* ------------------------------------------------------------------ */}
      {/* Header                                                               */}
      {/* ------------------------------------------------------------------ */}
      <PageHeader
        title="Cell-Type & Sex Resolution"
        description="Deconvolution attribution, SCENIC+ regulon analysis, and a sex-stratified pooled (cohort-adjusted) analysis reveal cell-type-specific and sex-dimorphic disease signatures."
      />

      {/* ------------------------------------------------------------------ */}
      {/* Section 1: Deconvolution Attribution                                 */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Deconvolution Attribution
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          MuSiC deconvolution reveals that 65.1% (975/1,497) of classifiable
          significant DEGs are hepatocyte-intrinsic. BayesPrism
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
                style={{ background: HEPATO_COLOR, opacity: 0.85 }}
              />
              <span>
                <span className="font-semibold">Hepatocyte-intrinsic</span>{" "}
                <span className="text-muted-foreground">
                  &mdash; 975 DEGs (65.1%)
                </span>
              </span>
            </div>
            <div className="flex items-center gap-2">
              <span
                className="inline-block h-3 w-3 rounded-sm"
                style={{ background: CONTROL, opacity: 0.55 }}
              />
              <span>
                <span className="font-semibold">Immune / stromal</span>{" "}
                <span className="text-muted-foreground">
                  &mdash; 522 DEGs (34.9%)
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
                        href={`#/gene?symbol=${encodeURIComponent(row.tf)}`}
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
                        style={{
                          color: isUp
                            ? "var(--color-effect-up)"
                            : "var(--color-effect-down)",
                        }}
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
          Sex-stratified interaction model (group_binary &times; inferred_sex),
          LVQW-fixed with dataset as a fixed effect (C2), reveals 8
          sex-dimorphic DEGs (interaction padj &lt; 0.05) &mdash; down from 290
          before fixing the dataset term. Direction concordance with the prior
          model is 0.80&ndash;0.84; only interaction significance collapses.
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
              <span
                className="inline-block h-2.5 w-2.5 rounded-full"
                style={{ backgroundColor: MALE_COLOR }}
              />
              Male-biased (7)
            </span>
            <span className="flex items-center gap-1.5">
              <span
                className="inline-block h-2.5 w-2.5 rounded-full"
                style={{ backgroundColor: FEMALE_COLOR }}
              />
              Female-biased (1)
            </span>
          </div>
        </div>

        {/* Methodology note */}
        <div className="rounded-lg border border-amber-300 bg-amber-50 px-5 py-4 dark:border-amber-700 dark:bg-amber-950/30">
          <p className="text-sm font-semibold text-amber-900 dark:text-amber-200">
            Methodology Note
          </p>
          <p className="mt-1 text-sm text-amber-800 dark:text-amber-300">
            An earlier sex-stratified model (dream, random dataset effect)
            reported 290 dimorphic genes with a female-enriched-progressor
            narrative. Fixing the dataset term to match the canonical C2 bulk
            model collapses this to 8 genes (7 male-biased / 1 female-biased).
            That earlier narrative is retired &mdash; at n&nbsp;=&nbsp;8, no
            COLOC-enrichment or -depletion claim can be made with statistical
            power; see the gene page for each of the 8 genes&apos; individual
            evidence.
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
    </PageContainer>
  );
}
