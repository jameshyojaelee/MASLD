"use client";

import { HashLink as Link } from "@/components/hash-link";
import { useHashLocation } from "@/lib/hash-router";
import { Search } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Separator } from "@/components/ui/separator";
import { ThemeToggle } from "@/components/theme-toggle";
import { useAppStore } from "@/lib/store";
import { NAV, isNavActive } from "@/lib/nav";
import { cn } from "@/lib/utils";

/** The MASLD wordmark + logo lockup. */
export function SidebarBrand() {
  return (
    <Link href="/" className="flex items-center gap-2">
      <div className="flex size-7 items-center justify-center rounded-md bg-primary text-primary-foreground text-xs font-bold">
        M
      </div>
      <span className="text-sm font-semibold tracking-tight">MASLD Gene Catalog</span>
    </Link>
  );
}

/** Command-palette trigger. */
export function SidebarSearchButton() {
  const { setCommandOpen } = useAppStore();
  return (
    <Button
      variant="outline"
      size="sm"
      className="w-full justify-start gap-2 text-muted-foreground"
      onClick={() => setCommandOpen(true)}
    >
      <Search className="size-3.5" aria-hidden />
      <span className="text-xs">Search genes...</span>
      <kbd className="ml-auto rounded border border-border bg-muted px-1 text-[10px] font-medium text-muted-foreground">
        {"⌘"}K
      </kbd>
    </Button>
  );
}

/**
 * The grouped navigation list, driven entirely by the NAV config. Shared by the
 * desktop sidebar and the mobile drawer. `onNavigate` lets the drawer close
 * itself when a link is chosen.
 */
export function SidebarNav({ onNavigate }: { onNavigate?: () => void }) {
  const { path: pathname } = useHashLocation();
  return (
    <nav className="flex-1 overflow-y-auto px-3 pb-3">
      {NAV.map((section) => (
        <div key={section.title} className="mb-3">
          <p className="mb-1 px-2 text-[10px] font-semibold uppercase tracking-wider text-muted-foreground">
            {section.title}
          </p>
          <ul className="space-y-0.5">
            {section.items.map((item) => {
              const active = isNavActive(pathname, item.href);
              const Icon = item.icon;
              return (
                <li key={item.href}>
                  <Link
                    href={item.href}
                    onClick={onNavigate}
                    className={cn(
                      "flex items-center gap-2 rounded-md px-2 py-1.5 text-sm transition-colors",
                      active
                        ? "bg-sidebar-accent text-sidebar-accent-foreground font-medium"
                        : "text-muted-foreground hover:bg-sidebar-accent hover:text-sidebar-accent-foreground"
                    )}
                  >
                    <Icon className="size-4 shrink-0" aria-hidden />
                    <span className="flex-1 truncate">{item.label}</span>
                    {item.flag && (
                      <span
                        title="Under review (keep / retire)"
                        aria-label="Under review"
                        className="size-1.5 shrink-0 rounded-full bg-muted-foreground/60"
                      />
                    )}
                  </Link>
                </li>
              );
            })}
          </ul>
        </div>
      ))}
    </nav>
  );
}

/** Desktop sidebar. Hidden below the `md` breakpoint (mobile uses the drawer). */
export function Sidebar() {
  return (
    <aside className="sticky top-0 hidden h-screen w-56 shrink-0 flex-col border-r border-sidebar-border bg-sidebar text-sidebar-foreground md:flex">
      <div className="px-4 py-4">
        <SidebarBrand />
      </div>

      <Separator />

      <div className="px-3 py-3">
        <SidebarSearchButton />
      </div>

      <Separator />

      <SidebarNav />

      <Separator />

      <div className="flex items-center justify-between px-4 py-3">
        <span className="text-[10px] font-medium text-muted-foreground">v2.0</span>
        <ThemeToggle />
      </div>
    </aside>
  );
}
