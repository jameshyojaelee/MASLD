"use client";

import { useState } from "react";
import { Menu } from "lucide-react";
import {
  Sidebar,
  SidebarBrand,
  SidebarNav,
  SidebarSearchButton,
} from "@/components/sidebar";
import { ThemeToggle } from "@/components/theme-toggle";
import { Button } from "@/components/ui/button";
import { Separator } from "@/components/ui/separator";
import { Sheet, SheetContent, SheetTitle } from "@/components/ui/sheet";

/**
 * Responsive application shell.
 *  - md+  : persistent left Sidebar + scrolling main.
 *  - < md : Sidebar collapses; a top bar exposes a Menu button that opens the
 *           navigation in a left Sheet drawer.
 */
export function AppShell({ children }: { children: React.ReactNode }) {
  const [open, setOpen] = useState(false);

  return (
    <div className="flex min-h-screen">
      <Sidebar />

      <div className="flex min-w-0 flex-1 flex-col">
        {/* Mobile top bar */}
        <header className="sticky top-0 z-40 flex items-center gap-2 border-b border-border bg-background/80 px-4 py-2 backdrop-blur md:hidden">
          <Button
            variant="ghost"
            size="icon-sm"
            onClick={() => setOpen(true)}
            aria-label="Open navigation"
          >
            <Menu className="size-5" aria-hidden />
          </Button>
          <SidebarBrand />
          <div className="ml-auto flex items-center gap-1">
            <ThemeToggle />
          </div>
        </header>

        <main className="flex-1 overflow-y-auto">{children}</main>
      </div>

      {/* Mobile navigation drawer */}
      <Sheet open={open} onOpenChange={setOpen}>
        <SheetContent
          side="left"
          className="flex w-72 flex-col gap-0 p-0 text-sidebar-foreground"
        >
          <SheetTitle className="sr-only">Navigation</SheetTitle>
          <div className="px-4 py-4">
            <SidebarBrand />
          </div>
          <Separator />
          <div className="px-3 py-3">
            <SidebarSearchButton />
          </div>
          <Separator />
          <div className="flex flex-1 flex-col overflow-y-auto">
            <SidebarNav onNavigate={() => setOpen(false)} />
          </div>
        </SheetContent>
      </Sheet>
    </div>
  );
}
