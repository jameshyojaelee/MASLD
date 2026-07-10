/**
 * Ambient "aurora" backdrop — decorative brand chrome for the landing hero.
 * Two-to-three large, soft, blurred radial-gradient blobs in the brand palette
 * (teal-blue + violet, with a warm amber counterpoint) that drift slowly behind
 * the constellation and headline, reading as depth rather than motion.
 *
 * CSS-only: each blob runs an `aurora-drift-*` keyframe (globals.css) with its
 * own duration / negative delay / start position so they wander independently,
 * and freezes to a static wash under `prefers-reduced-motion`. Dark-first, low
 * opacity, pointer-events-none, aria-hidden. Encodes no data — no scale, no
 * legend; the colors are chrome, never a data mark.
 */

type Blob = {
  /** Radial-gradient core color — decorative chrome, never a data mark. */
  color: string;
  /** Position + size utilities (blobs sit off-center so drift stays diffuse). */
  className: string;
  opacity: number;
  animationName: string;
  animationDuration: string;
  /** Negative delay desyncs the blobs from a shared frame-0 start. */
  animationDelay: string;
};

const BLOBS: readonly Blob[] = [
  {
    color: "var(--primary)",
    className: "left-[-10%] top-[-15%] h-[40rem] w-[40rem]",
    opacity: 0.24,
    animationName: "aurora-drift-a",
    animationDuration: "22s",
    animationDelay: "0s",
  },
  {
    color: "var(--accent-brand)",
    className: "right-[-12%] top-[6%] h-[34rem] w-[34rem]",
    opacity: 0.22,
    animationName: "aurora-drift-b",
    animationDuration: "28s",
    animationDelay: "-6s",
  },
  {
    // Warm amber counterpoint to the cool brand duo (fixed hue, both themes).
    color: "oklch(0.80 0.13 68)",
    className: "bottom-[-18%] left-[28%] h-[30rem] w-[30rem]",
    opacity: 0.16,
    animationName: "aurora-drift-c",
    animationDuration: "34s",
    animationDelay: "-12s",
  },
];

export function Aurora() {
  return (
    <div
      aria-hidden
      className="pointer-events-none absolute inset-0 -z-[20] overflow-hidden"
    >
      {BLOBS.map((b, i) => (
        <div
          key={i}
          className={`aurora-blob ${b.className}`}
          style={{
            background: `radial-gradient(circle at center, ${b.color} 0%, transparent 70%)`,
            opacity: b.opacity,
            animationName: b.animationName,
            animationDuration: b.animationDuration,
            animationDelay: b.animationDelay,
          }}
        />
      ))}
    </div>
  );
}

export default Aurora;
