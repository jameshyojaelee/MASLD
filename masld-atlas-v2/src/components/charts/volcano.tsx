"use client";

/**
 * Volcano plot: log2 fold-change (x) vs −log10(FDR) (y).
 *
 * Direction is encoded by hue (blue = down, red = up) via the diverging scale,
 * and significance is ALSO encoded by opacity — non-significant points fall back
 * to control gray at reduced opacity (never color alone). Vertical guides mark
 * the ±logFC effect-size floor; a horizontal guide marks the FDR threshold.
 *
 * Interactivity (all opt-in, non-breaking):
 *   - `enableBrush`   drag a rectangle to zoom; double-click resets.
 *   - `autoLabelTopN` auto-annotate the N largest-effect significant genes.
 *   - click a point   pins/unpins its label (unless `onPointClick` is set, in
 *                      which case a click navigates instead).
 *   - staggered mark mount (size-gated) — animates element opacity only, never
 *     the resting significance-opacity encoding.
 */

import { useCallback, useId, useMemo, useRef, useState } from "react";
import { scaleLinear } from "@visx/scale";
import { Line } from "@visx/shape";
import { Text } from "@visx/text";
import { m, useReducedMotion } from "framer-motion";
import { getColor, CONTROL, NONSIG_OPACITY } from "@/lib/palette";
import { DEG_LFC, DEG_FDR } from "@/lib/atlas-constants";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import {
  BottomAxis,
  LeftAxis,
  GridRows,
  GridCols,
  AXIS_LINE,
  MARK_STAGGER_MAX,
  markItemVariants,
  MarkGroup,
} from "./chart-primitives";

export interface VolcanoPoint {
  symbol: string;
  logFC: number;
  /** Adjusted p-value / FDR; internally transformed to −log10. */
  padj: number;
  /** Explicit significance; if omitted it is derived from the thresholds. */
  sig?: boolean;
}

export interface VolcanoProps {
  data: VolcanoPoint[];
  /** Effect-size floor for the vertical guides (default = canonical DEG LFC). */
  lfcThreshold?: number;
  /** FDR threshold for the horizontal guide + derived significance. */
  fdrThreshold?: number;
  /** Symbols to annotate with a text label. */
  highlightSymbols?: string[];
  /** Auto-annotate the N largest-effect significant genes (in addition to
   * `highlightSymbols` and any pinned points). Default 0 (off). */
  autoLabelTopN?: number;
  /** Enable drag-rectangle zoom (double-click to reset). Default off. */
  enableBrush?: boolean;
  title?: React.ReactNode;
  caption?: React.ReactNode;
  height?: number;
  ariaLabel?: string;
  onPointClick?: (symbol: string) => void;
}

const MIN_LOG10 = -Math.log10(1); // 0
const EPS = 1e-300;

function negLog10(p: number): number {
  return -Math.log10(Math.max(p, EPS));
}

type Domain = [number, number];
interface ZoomState {
  x: Domain;
  y: Domain;
}

export function Volcano({
  data,
  lfcThreshold = DEG_LFC,
  fdrThreshold = DEG_FDR,
  highlightSymbols,
  autoLabelTopN = 0,
  enableBrush = false,
  title,
  caption,
  height = 360,
  ariaLabel,
  onPointClick,
}: VolcanoProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<VolcanoPoint>();
  const clipId = useId().replace(/:/g, "");
  const reduce = useReducedMotion();

  const [pinned, setPinned] = useState<Set<string>>(new Set());
  const [zoom, setZoom] = useState<ZoomState | null>(null);
  // Live brush rectangle in plot-local px (null when not brushing).
  const [brush, setBrush] = useState<{ x0: number; y0: number; x1: number; y1: number } | null>(
    null
  );
  const brushStart = useRef<{ x: number; y: number } | null>(null);
  const captureRectRef = useRef<SVGRectElement>(null);

  const togglePin = useCallback((symbol: string) => {
    setPinned((prev) => {
      const next = new Set(prev);
      if (next.has(symbol)) next.delete(symbol);
      else next.add(symbol);
      return next;
    });
  }, []);

  const points = useMemo(
    () =>
      data.map((d) => {
        const y = negLog10(d.padj);
        const sig =
          d.sig ?? (d.padj < fdrThreshold && Math.abs(d.logFC) >= lfcThreshold);
        return { ...d, y, sig };
      }),
    [data, fdrThreshold, lfcThreshold]
  );

  // Labels = explicit highlights ∪ pinned ∪ top-N largest-effect significant.
  const labelSet = useMemo(() => {
    const s = new Set(highlightSymbols ?? []);
    for (const p of pinned) s.add(p);
    if (autoLabelTopN > 0) {
      [...points]
        .filter((p) => p.sig)
        .sort((a, b) => Math.abs(b.logFC) - Math.abs(a.logFC))
        .slice(0, autoLabelTopN)
        .forEach((p) => s.add(p.symbol));
    }
    return s;
  }, [highlightSymbols, pinned, autoLabelTopN, points]);

  const fullXDomain = useMemo<Domain>(() => {
    const maxAbs = Math.max(1, ...points.map((p) => Math.abs(p.logFC)));
    const pad = maxAbs * 0.05;
    return [-(maxAbs + pad), maxAbs + pad];
  }, [points]);

  const fullYMax = useMemo(
    () => Math.max(negLog10(fdrThreshold) * 1.2, ...points.map((p) => p.y), 1),
    [points, fdrThreshold]
  );

  const xDomain = zoom?.x ?? fullXDomain;
  const yDomain: Domain = zoom?.y ?? [MIN_LOG10, fullYMax];
  const zoomed = zoom != null;

  // Gate only on mark count + reduced motion (NOT zoom) so the mark container
  // doesn't switch element types on zoom and re-trigger the mount reveal.
  const staggerMarks = !reduce && points.length <= MARK_STAGGER_MAX;

  return (
    <div ref={wrapperRef} className="relative w-full">
      {enableBrush && zoomed && (
        <button
          type="button"
          onClick={() => setZoom(null)}
          className="absolute right-2 top-1 z-20 rounded border border-border bg-background/90 px-2 py-0.5 text-[11px] text-muted-foreground shadow-sm backdrop-blur hover:text-foreground"
        >
          Reset zoom
        </button>
      )}
      <ChartFrame
        title={title}
        caption={caption}
        height={height}
        ariaLabel={ariaLabel ?? "Volcano plot"}
        margin={{ top: 12, right: 16, bottom: 44, left: 52 }}
      >
        {({ innerWidth, innerHeight }) => {
          const xScale = scaleLinear({ domain: xDomain, range: [0, innerWidth] });
          const yScale = scaleLinear({
            domain: yDomain,
            range: [innerHeight, 0],
            nice: !zoomed,
          });
          const xTicks = xScale.ticks(7);
          const yTicks = yScale.ticks(5);
          const fdrY = yScale(negLog10(fdrThreshold));

          // Brush pointer handlers (capture local scales via this closure).
          const localXY = (e: React.PointerEvent) => {
            const rect = captureRectRef.current?.getBoundingClientRect();
            if (!rect) return null;
            return {
              x: Math.max(0, Math.min(innerWidth, e.clientX - rect.left)),
              y: Math.max(0, Math.min(innerHeight, e.clientY - rect.top)),
            };
          };
          const onBrushDown = (e: React.PointerEvent) => {
            const p = localXY(e);
            if (!p) return;
            brushStart.current = p;
            setBrush({ x0: p.x, y0: p.y, x1: p.x, y1: p.y });
            captureRectRef.current?.setPointerCapture(e.pointerId);
          };
          const onBrushMove = (e: React.PointerEvent) => {
            if (!brushStart.current) return;
            const p = localXY(e);
            if (!p) return;
            setBrush({ x0: brushStart.current.x, y0: brushStart.current.y, x1: p.x, y1: p.y });
          };
          const onBrushUp = (e: React.PointerEvent) => {
            const start = brushStart.current;
            brushStart.current = null;
            captureRectRef.current?.releasePointerCapture?.(e.pointerId);
            const p = localXY(e);
            setBrush(null);
            if (!start || !p) return;
            const dx = Math.abs(p.x - start.x);
            const dy = Math.abs(p.y - start.y);
            if (dx < 6 || dy < 6) return; // ignore stray clicks
            const nx: Domain = [
              xScale.invert(Math.min(start.x, p.x)),
              xScale.invert(Math.max(start.x, p.x)),
            ];
            // y range is inverted (top = high value)
            const ny: Domain = [
              yScale.invert(Math.max(start.y, p.y)),
              yScale.invert(Math.min(start.y, p.y)),
            ];
            setZoom({ x: nx, y: ny });
          };

          const marks = points.map((p, i) => {
            const cx = xScale(p.logFC);
            const cy = yScale(p.y);
            if (zoomed && (cx < 0 || cx > innerWidth || cy < 0 || cy > innerHeight)) {
              return null;
            }
            const { fill, opacity } = getColor(p.logFC, p.sig);
            const isHi = labelSet.has(p.symbol);
            const common = {
              cx,
              cy,
              r: isHi ? 4 : 2.6,
              fill: p.sig ? fill : CONTROL,
              fillOpacity: p.sig ? opacity : NONSIG_OPACITY,
              stroke: isHi ? "var(--color-foreground)" : "none",
              strokeWidth: isHi ? 1 : 0,
              className: "cursor-pointer",
              onMouseMove: (e: React.MouseEvent) => show(e, p),
              onMouseLeave: hide,
              onClick: () => (onPointClick ? onPointClick(p.symbol) : togglePin(p.symbol)),
            };
            return staggerMarks ? (
              <m.circle
                key={`${p.symbol}-${i}`}
                custom={i}
                variants={markItemVariants}
                initial="hidden"
                animate="show"
                {...common}
              />
            ) : (
              <circle key={`${p.symbol}-${i}`} {...common} />
            );
          });

          return (
            <>
              <defs>
                <clipPath id={`clip-${clipId}`}>
                  <rect x={0} y={0} width={innerWidth} height={innerHeight} />
                </clipPath>
              </defs>

              <GridRows scale={yScale} ticks={yTicks} width={innerWidth} />
              <GridCols scale={xScale} ticks={xTicks} height={innerHeight} />

              {/* Effect-size floor guides */}
              {[-lfcThreshold, lfcThreshold].map((t) => (
                <Line
                  key={t}
                  from={{ x: xScale(t), y: 0 }}
                  to={{ x: xScale(t), y: innerHeight }}
                  stroke={AXIS_LINE}
                  strokeDasharray="3 3"
                  strokeWidth={1}
                />
              ))}
              {/* FDR guide */}
              <Line
                from={{ x: 0, y: fdrY }}
                to={{ x: innerWidth, y: fdrY }}
                stroke={AXIS_LINE}
                strokeDasharray="3 3"
                strokeWidth={1}
              />

              {/* Brush capture surface (below the marks so points keep hover) */}
              {enableBrush && (
                <rect
                  ref={captureRectRef}
                  x={0}
                  y={0}
                  width={innerWidth}
                  height={innerHeight}
                  fill="transparent"
                  style={{ cursor: "crosshair", touchAction: "none" }}
                  onPointerDown={onBrushDown}
                  onPointerMove={onBrushMove}
                  onPointerUp={onBrushUp}
                  onDoubleClick={() => setZoom(null)}
                />
              )}

              {/* Points */}
              <g clipPath={`url(#clip-${clipId})`}>
                {staggerMarks ? marks : <MarkGroup>{marks}</MarkGroup>}
              </g>

              {/* Highlight labels */}
              <g clipPath={`url(#clip-${clipId})`}>
                {points
                  .filter((p) => labelSet.has(p.symbol))
                  .map((p, i) => {
                    const cx = xScale(p.logFC);
                    const cy = yScale(p.y);
                    if (zoomed && (cx < 0 || cx > innerWidth || cy < 0 || cy > innerHeight))
                      return null;
                    return (
                      <Text
                        key={`lbl-${p.symbol}-${i}`}
                        x={cx}
                        y={cy - 7}
                        fontSize={10}
                        textAnchor="middle"
                        fill="var(--color-foreground)"
                        fontStyle="italic"
                      >
                        {p.symbol}
                      </Text>
                    );
                  })}
              </g>

              {/* Live brush rectangle (chrome; never steals pointer events) */}
              {brush && (
                <rect
                  x={Math.min(brush.x0, brush.x1)}
                  y={Math.min(brush.y0, brush.y1)}
                  width={Math.abs(brush.x1 - brush.x0)}
                  height={Math.abs(brush.y1 - brush.y0)}
                  fill={CONTROL}
                  fillOpacity={0.12}
                  stroke={AXIS_LINE}
                  strokeDasharray="3 3"
                  pointerEvents="none"
                />
              )}

              <BottomAxis
                scale={xScale}
                top={innerHeight}
                numTicks={7}
                label="log₂ fold-change"
              />
              <LeftAxis scale={yScale} numTicks={5} label="−log₁₀ FDR" />
            </>
          );
        }}
      </ChartFrame>

      {tooltip && (
        <ChartTooltip left={tooltip.left} top={tooltip.top}>
          <div className="font-medium italic">{tooltip.data.symbol}</div>
          <div>logFC {tooltip.data.logFC.toFixed(2)}</div>
          <div>FDR {tooltip.data.padj.toExponential(1)}</div>
        </ChartTooltip>
      )}
    </div>
  );
}
