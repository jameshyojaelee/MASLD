"use client";

/**
 * Tiny element-size hook shared by the landing mini-viz previews.
 *
 * Previews draw into a measured pixel box (rather than a scaled SVG viewBox) so
 * point radii, canvas dots, and axis text keep a constant physical size no
 * matter how large the containing bento cell is. Returns `[ref, {width,
 * height}]`; size is `0` until the first ResizeObserver callback, so callers
 * should render nothing (or a skeleton) while width is 0.
 */

import { useEffect, useRef, useState } from "react";

export interface Size {
  width: number;
  height: number;
}

export function useMeasure<T extends HTMLElement = HTMLDivElement>() {
  const ref = useRef<T>(null);
  const [size, setSize] = useState<Size>({ width: 0, height: 0 });

  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    // Seed synchronously so the first paint after mount already has a box.
    setSize({ width: el.clientWidth, height: el.clientHeight });
    if (typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver((entries) => {
      const cr = entries[0]?.contentRect;
      if (cr) setSize({ width: cr.width, height: cr.height });
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  return [ref, size] as const;
}
