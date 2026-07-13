import type { NextConfig } from "next";

// NOTE on routing (Hugging Face Static Space):
// HF's static server resolves an extensionless request path P by appending
// `.html` (`/single-cell` -> `single-cell.html`), but it does NOT serve
// `<dir>/index.html` for a bare `<dir>/` — a trailing-slash path apex-redirects
// to huggingface.co/<dir>. So `trailingSlash: true` (which emits
// `out/<route>/index.html`) breaks every sub-route hard-refresh / deep-link.
// Leaving trailingSlash at its default (false) makes `output: export` emit
// `out/<route>.html`, which HF serves at `/<route>` via `.html`-append.
// Consequence: all internal gene links use `/gene?symbol=…` (NO trailing slash).
const nextConfig: NextConfig = {
  output: "export",
  images: {
    unoptimized: true,
  },
};

export default nextConfig;
