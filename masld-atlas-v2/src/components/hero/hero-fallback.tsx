/**
 * Static SVG constellation — the hero's zero-JS first paint and its permanent
 * stand-in under reduced motion, no canvas, or Save-Data. Fixed dots + faint
 * links, painted with the brand gradient (decorative chrome). No graph library,
 * no interaction. Positions are hand-placed (deterministic) so it reads as a
 * constellation without implying any quantitative relationship.
 */

// [x, y, r] in the 100 x 60 viewBox.
const STARS: readonly [number, number, number][] = [
  [8, 12, 1.4], [16, 28, 2.6], [12, 46, 1.2], [24, 16, 1.8], [30, 38, 2.2],
  [22, 52, 1.5], [38, 10, 1.3], [42, 30, 3.0], [36, 50, 1.6], [50, 20, 1.9],
  [54, 42, 2.4], [48, 54, 1.2], [62, 14, 1.7], [66, 34, 2.8], [60, 50, 1.4],
  [74, 24, 2.0], [78, 44, 1.6], [72, 56, 1.2], [86, 16, 1.5], [90, 36, 2.3],
  [84, 52, 1.3], [94, 26, 1.6], [20, 38, 1.1], [58, 28, 1.2],
];

const LINKS: readonly [number, number][] = [
  [0, 1], [1, 3], [3, 6], [1, 4], [4, 7], [4, 22], [7, 9], [9, 12], [12, 15],
  [15, 18], [18, 21], [7, 13], [13, 16], [16, 19], [10, 13], [10, 23], [4, 10],
  [1, 2], [2, 5], [8, 11], [14, 17], [17, 20], [19, 20], [6, 9],
];

export function HeroFallback() {
  return (
    <svg
      aria-hidden
      className="h-full w-full"
      viewBox="0 0 100 60"
      preserveAspectRatio="xMidYMid slice"
    >
      <defs>
        <linearGradient id="hero-constellation-grad" x1="0" y1="0" x2="1" y2="1">
          <stop offset="0" stopColor="var(--primary)" />
          <stop offset="1" stopColor="var(--accent-brand)" />
        </linearGradient>
      </defs>
      <g stroke="url(#hero-constellation-grad)" strokeWidth="0.15" opacity="0.35">
        {LINKS.map(([a, b], i) => (
          <line
            key={i}
            x1={STARS[a][0]}
            y1={STARS[a][1]}
            x2={STARS[b][0]}
            y2={STARS[b][1]}
          />
        ))}
      </g>
      <g fill="url(#hero-constellation-grad)">
        {STARS.map(([x, y, r], i) => (
          <circle key={`g${i}`} cx={x} cy={y} r={r * 2.6} opacity="0.12" />
        ))}
        {STARS.map(([x, y, r], i) => (
          <circle key={`s${i}`} cx={x} cy={y} r={r} opacity="0.85" />
        ))}
      </g>
    </svg>
  );
}

export default HeroFallback;
