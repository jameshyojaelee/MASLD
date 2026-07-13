"use client";

/**
 * App-wide motion provider.
 *
 * `LazyMotion` + `domAnimation` loads only the DOM animation feature bundle
 * (~15 KB) instead of the full `motion.*` runtime — all animated components
 * import the lightweight `m` component (never `motion`) so this tree-shakes.
 *
 * `MotionConfig reducedMotion="user"` makes EVERY framer animation honor the
 * user's `prefers-reduced-motion` setting globally, so individual components
 * don't each have to gate their variants.
 */

import { LazyMotion, domAnimation, MotionConfig } from "framer-motion";
import type { ReactNode } from "react";

export function MotionProvider({ children }: { children: ReactNode }) {
  return (
    <LazyMotion features={domAnimation} strict={false}>
      <MotionConfig reducedMotion="user">{children}</MotionConfig>
    </LazyMotion>
  );
}
