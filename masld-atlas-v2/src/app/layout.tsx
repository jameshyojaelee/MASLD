import type { Metadata } from "next";
import { Inter, JetBrains_Mono, Space_Grotesk } from "next/font/google";
import "./globals.css";
import { AppShell } from "@/components/app-shell";
import { MotionProvider } from "@/components/motion/motion-provider";
import { HashRouterProvider } from "@/lib/hash-router";
import { RouteSwitch } from "@/components/route-switch";
import { GeneSearch } from "@/components/gene-search";
import { GeneComparePanel } from "@/components/gene-compare-panel";

const inter = Inter({
  subsets: ["latin"],
  display: "swap",
  variable: "--font-inter",
});

const jetbrainsMono = JetBrains_Mono({
  subsets: ["latin"],
  display: "swap",
  variable: "--font-jetbrains-mono",
});

// Display face for headlines + hero numerals (variable font; self-hosted).
const spaceGrotesk = Space_Grotesk({
  subsets: ["latin"],
  display: "swap",
  variable: "--font-space-grotesk",
});

export const metadata: Metadata = {
  title: "MASLD Gene Catalog — legacy pre-Resource preview",
  description:
    "Legacy pre-Resource implementation; publication deployment is closed pending the MASLD Gene Catalog release",
};

/**
 * Pre-hydration theme script. Runs synchronously as the first child of <body>,
 * before body content paints, so the correct theme is applied with no flash.
 * Reads the persisted choice (localStorage). The app is dark-first, so it
 * defaults to dark unless the user has explicitly chosen light.
 * zustand remains the runtime source of truth once React hydrates.
 */
const THEME_INIT = `
(function () {
  try {
    var s = localStorage.getItem('masld-atlas-theme');
    var dark = s !== 'light'; // dark by default; only an explicit 'light' opts out
    document.documentElement.classList.toggle('dark', dark);
  } catch (e) {}
})();
`;

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html
      lang="en"
      className={`${inter.variable} ${jetbrainsMono.variable} ${spaceGrotesk.variable}`}
      suppressHydrationWarning
    >
      <body className="min-h-screen antialiased" suppressHydrationWarning>
        <script dangerouslySetInnerHTML={{ __html: THEME_INIT }} />
        <MotionProvider>
          <HashRouterProvider>
            <AppShell>
              <div
                role="alert"
                className="border-b border-amber-500/50 bg-amber-500/10 px-4 py-2 text-center text-xs text-amber-200"
              >
                Legacy pre-Resource preview — not publication-ready. Causal-gene,
                longitudinal-progression, universal-ranking, and Cas13-target
                views are retired; the publication portal will default to
                the MASLD Gene Catalog and gene-specific follow-up guidance.
              </div>
              <RouteSwitch>{children}</RouteSwitch>
            </AppShell>
            <GeneComparePanel />
            <GeneSearch />
          </HashRouterProvider>
        </MotionProvider>
      </body>
    </html>
  );
}
