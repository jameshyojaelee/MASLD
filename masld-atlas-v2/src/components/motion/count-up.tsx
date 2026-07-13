"use client";

/**
 * Animated count-up for headline metrics.
 *
 * Uses framer's imperative `animate()` to tween a raw number. Under reduced
 * motion (or when `animateOnMount` is false) it renders the final formatted
 * value immediately — no tween, no layout jitter. Pass ONLY real numbers here;
 * string/ReactNode values should render directly without this component.
 */

import { useEffect, useRef, useState } from "react";
import { animate, useReducedMotion } from "framer-motion";
import { MOTION } from "@/lib/motion";

export function CountUp({
  value,
  format = (n: number) => Math.round(n).toLocaleString(),
  durationMs = 900,
  className,
  animateOnMount = true,
}: {
  value: number;
  format?: (n: number) => string;
  durationMs?: number;
  className?: string;
  animateOnMount?: boolean;
}) {
  const reduce = useReducedMotion();
  // Initial render MUST NOT depend on `reduce`: it is false during SSR/prerender
  // but can be true on the client's first paint under reduced motion, which
  // desyncs the server and client text and triggers a hydration mismatch
  // (React #418). Start deterministically at 0 (animated) / value (static), then
  // let the effect jump-to-value (reduced motion) or animate.
  const [display, setDisplay] = useState<number>(animateOnMount ? 0 : value);
  const fromRef = useRef<number>(animateOnMount ? 0 : value);

  useEffect(() => {
    if (reduce || !animateOnMount) {
      setDisplay(value);
      fromRef.current = value;
      return;
    }
    const controls = animate(fromRef.current, value, {
      duration: durationMs / 1000,
      ease: MOTION.easeOutExpo,
      onUpdate: (v) => setDisplay(v),
    });
    fromRef.current = value;
    return () => controls.stop();
  }, [value, reduce, durationMs, animateOnMount]);

  return <span className={className}>{format(display)}</span>;
}
