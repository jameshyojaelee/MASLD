/**
 * <MarqueeRow> — dependency-free, seamless CSS marquee.
 *
 * Renders its children TWICE inside an overflow-hidden viewport; the track
 * slides via a pure CSS `transform` keyframe (`.animate-marquee` in globals.css)
 * from 0 to -50%, i.e. exactly one copy-width, so the loop is gapless. No
 * framer / JS animation is involved (the marquee doctrine forbids framer
 * `layout`; CSS transform is the sanctioned path).
 *
 * The second copy is `aria-hidden` + `inert` so screen readers and tab order
 * see the content only once. Under `prefers-reduced-motion: reduce` the CSS
 * utility freezes the track (no translate) and the overflow clip leaves a single
 * static row visible. `pauseOnHover` and the caller-driven `paused` prop both
 * map onto `animation-play-state: paused`.
 */

import type { CSSProperties, ReactNode } from "react";
import { cn } from "@/lib/utils";

export interface MarqueeRowProps {
  children: ReactNode;
  /** Seconds for one full loop. Larger = slower. Default ~40s (slow). */
  durationSec?: number;
  /** Pause the scroll while the pointer is over the row. Default true. */
  pauseOnHover?: boolean;
  /** Externally freeze the scroll (e.g. when offscreen). Default false. */
  paused?: boolean;
  className?: string;
}

export function MarqueeRow({
  children,
  durationSec = 40,
  pauseOnHover = true,
  paused = false,
  className,
}: MarqueeRowProps) {
  return (
    <div
      className={cn(
        "overflow-hidden",
        pauseOnHover && "marquee-hoverpause",
        className
      )}
    >
      <div
        className="animate-marquee flex w-max flex-nowrap"
        data-paused={paused ? "true" : "false"}
        style={{ "--marquee-duration": `${durationSec}s` } as CSSProperties}
      >
        <div className="flex shrink-0 flex-nowrap items-center">{children}</div>
        <div
          className="flex shrink-0 flex-nowrap items-center"
          aria-hidden
          // React 19 boolean `inert`: drops the duplicate copy from a11y + tab order.
          inert
        >
          {children}
        </div>
      </div>
    </div>
  );
}

export default MarqueeRow;
