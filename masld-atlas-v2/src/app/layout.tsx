import type { Metadata } from "next";
import { Inter } from "next/font/google";
import "./globals.css";
import { Sidebar } from "@/components/sidebar";
import { GeneSearch } from "@/components/gene-search";
import { GeneComparePanel } from "@/components/gene-compare-panel";

const inter = Inter({
  subsets: ["latin"],
  display: "swap",
  variable: "--font-inter",
});

export const metadata: Metadata = {
  title: "MASLD Atlas",
  description: "Multi-evidence transcriptomic atlas for metabolic dysfunction-associated steatotic liver disease",
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en" className={inter.variable} suppressHydrationWarning>
      <body className="min-h-screen antialiased">
        <div className="flex">
          <Sidebar />
          <main className="flex-1 overflow-y-auto">{children}</main>
        </div>
        <GeneComparePanel />
        <GeneSearch />
      </body>
    </html>
  );
}
