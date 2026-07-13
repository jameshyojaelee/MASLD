"use client";

/**
 * Dot plot: gene (row) × group (col), dot AREA encodes a fraction (e.g. % of
 * cells expressing) and dot COLOR encodes a mean (e.g. scaled mean expression)
 * on the sequential ramp. Both channels get a legend. The standard scRNA
 * marker-panel idiom.
 */

import { useMemo, useState } from "react";
import { scaleBand, scaleSqrt } from "@visx/scale";
import { Text } from "@visx/text";
import { sequentialColor } from "@/lib/palette";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import { SEQUENTIAL_RAMP } from "@/lib/palette";
import { AXIS_TEXT, MarkGroup } from "./chart-primitives";

export interface DotPlotPoint {
  row: string;
  col: string;
  /** Size channel, a fraction in [0, 1]. */
  size: number;
  /** Color channel, e.g. mean expression (mapped over its own domain). */
  color: number;
}

export interface DotPlotProps {
  data: DotPlotPoint[];
  rows?: string[];
  cols?: string[];
  /** Color domain; defaults to data min/max of `color`. */
  colorDomain?: [number, number];
  colorLabel?: string;
  sizeLabel?: string;
  valueFormat?: (v: number) => string;
  title?: React.ReactNode;
  caption?: React.ReactNode;
  height?: number;
  ariaLabel?: string;
}

function uniqueInOrder(values: string[]): string[] {
  const seen = new Set<string>();
  const out: string[] = [];
  for (const v of values) {
    if (!seen.has(v)) {
      seen.add(v);
      out.push(v);
    }
  }
  return out;
}

const KEY = (r: string, c: string) => `${r} ${c}`;

export function DotPlot({
  data,
  rows: rowsProp,
  cols: colsProp,
  colorDomain,
  colorLabel = "mean",
  sizeLabel = "fraction",
  valueFormat = (v) => v.toFixed(2),
  title,
  caption,
  height,
  ariaLabel,
}: DotPlotProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<DotPlotPoint>();
  const [hoverKey, setHoverKey] = useState<{ row: string; col: string } | null>(null);

  const rows = useMemo(
    () => rowsProp ?? uniqueInOrder(data.map((d) => d.row)),
    [rowsProp, data]
  );
  const cols = useMemo(
    () => colsProp ?? uniqueInOrder(data.map((d) => d.col)),
    [colsProp, data]
  );

  const [cMin, cMax] = useMemo(() => {
    if (colorDomain) return colorDomain;
    const vals = data.map((d) => d.color);
    return [Math.min(...vals), Math.max(...vals)] as [number, number];
  }, [colorDomain, data]);

  const colorFor = useMemo(() => {
    const span = cMax - cMin || 1;
    return (v: number) => sequentialColor((v - cMin) / span);
  }, [cMin, cMax]);

  const resolvedHeight =
    height ?? Math.min(560, Math.max(160, rows.length * 30 + 96));

  const rowLabelW = 96;
  const colLabelH = 64;

  return (
    <div ref={wrapperRef} className="relative w-full">
      <ChartFrame
        title={title}
        caption={caption}
        height={resolvedHeight}
        ariaLabel={ariaLabel ?? "Dot plot"}
        margin={{ top: colLabelH, right: 16, bottom: 8, left: rowLabelW }}
      >
        {({ innerWidth, innerHeight }) => {
          const xScale = scaleBand({ domain: cols, range: [0, innerWidth], padding: 0.1 });
          const yScale = scaleBand({ domain: rows, range: [0, innerHeight], padding: 0.1 });
          const maxR = Math.min(xScale.bandwidth(), yScale.bandwidth()) / 2 - 1;
          const rScale = scaleSqrt({ domain: [0, 1], range: [1, Math.max(3, maxR)] });

          return (
            <>
              <MarkGroup>
                {rows.map((r) =>
                  cols.map((c) => {
                    const d = data.find((p) => p.row === r && p.col === c);
                    if (!d) return null;
                    const cx = (xScale(c) ?? 0) + xScale.bandwidth() / 2;
                    const cy = (yScale(r) ?? 0) + yScale.bandwidth() / 2;
                    const dim =
                      hoverKey != null && r !== hoverKey.row && c !== hoverKey.col;
                    return (
                      <circle
                        key={KEY(r, c)}
                        cx={cx}
                        cy={cy}
                        r={rScale(Math.max(0, Math.min(1, d.size)))}
                        fill={colorFor(d.color)}
                        fillOpacity={dim ? 0.28 : 1}
                        stroke="var(--color-border)"
                        strokeWidth={0.5}
                        onMouseMove={(e) => {
                          show(e, d);
                          setHoverKey({ row: r, col: c });
                        }}
                        onMouseLeave={() => {
                          hide();
                          setHoverKey(null);
                        }}
                      />
                    );
                  })
                )}
              </MarkGroup>

              {rows.map((r) => (
                <Text
                  key={`r-${r}`}
                  x={-6}
                  y={(yScale(r) ?? 0) + yScale.bandwidth() / 2}
                  verticalAnchor="middle"
                  textAnchor="end"
                  fontSize={10}
                  fontStyle="italic"
                  fill={AXIS_TEXT}
                  fontWeight={hoverKey?.row === r ? 700 : 400}
                  width={rowLabelW - 10}
                >
                  {r}
                </Text>
              ))}

              {cols.map((c) => (
                <Text
                  key={`c-${c}`}
                  x={(xScale(c) ?? 0) + xScale.bandwidth() / 2}
                  y={-6}
                  verticalAnchor="middle"
                  textAnchor="start"
                  angle={-45}
                  fontSize={10}
                  fill={AXIS_TEXT}
                  fontWeight={hoverKey?.col === c ? 700 : 400}
                >
                  {c}
                </Text>
              ))}
            </>
          );
        }}
      </ChartFrame>

      {/* Combined size + color legend */}
      <div className="mt-2 flex flex-wrap items-center gap-x-6 gap-y-2 text-xs text-muted-foreground">
        <div className="flex items-center gap-2">
          <span>{sizeLabel}</span>
          {[0.1, 0.5, 1].map((f) => (
            <span key={f} className="inline-flex items-center gap-1">
              <svg width={16} height={16} aria-hidden>
                <circle
                  cx={8}
                  cy={8}
                  r={Math.max(2, f * 7)}
                  fill="var(--color-muted-foreground)"
                  fillOpacity={0.6}
                />
              </svg>
              <span className="tabular-nums">{f}</span>
            </span>
          ))}
        </div>
        <div className="flex items-center gap-2">
          <span>{colorLabel}</span>
          <span className="tabular-nums">{valueFormat(cMin)}</span>
          <span
            className="inline-block h-2.5 w-16 rounded-full"
            style={{
              background: `linear-gradient(to right, ${SEQUENTIAL_RAMP.join(", ")})`,
            }}
          />
          <span className="tabular-nums">{valueFormat(cMax)}</span>
        </div>
      </div>

      {tooltip && (
        <ChartTooltip left={tooltip.left} top={tooltip.top}>
          <div className="font-medium italic">{tooltip.data.row}</div>
          <div>{tooltip.data.col}</div>
          <div>
            {sizeLabel} {(tooltip.data.size * 100).toFixed(0)}%
          </div>
          <div>
            {colorLabel} {valueFormat(tooltip.data.color)}
          </div>
        </ChartTooltip>
      )}
    </div>
  );
}
