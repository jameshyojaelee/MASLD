"use client";

import { useCallback, useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import {
  Command,
  CommandEmpty,
  CommandGroup,
  CommandInput,
  CommandItem,
  CommandList,
} from "@/components/ui/command";
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
} from "@/components/ui/dialog";
import { Badge } from "@/components/ui/badge";
import { useAppStore } from "@/lib/store";
import { searchGenes } from "@/lib/search-index";
import type { GeneIndexEntry } from "@/lib/types";

export function GeneSearch() {
  const router = useRouter();
  const { commandOpen, setCommandOpen } = useAppStore();
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<GeneIndexEntry[]>([]);

  useEffect(() => {
    const down = (e: KeyboardEvent) => {
      if (e.key === "k" && (e.metaKey || e.ctrlKey)) {
        e.preventDefault();
        setCommandOpen(!commandOpen);
      }
    };
    document.addEventListener("keydown", down);
    return () => document.removeEventListener("keydown", down);
  }, [commandOpen, setCommandOpen]);

  useEffect(() => {
    if (!query || query.length < 2) {
      setResults([]);
      return;
    }
    let cancelled = false;
    searchGenes(query, 10).then((r) => {
      if (!cancelled) setResults(r);
    });
    return () => {
      cancelled = true;
    };
  }, [query]);

  const handleSelect = useCallback(
    (symbol: string) => {
      setCommandOpen(false);
      setQuery("");
      router.push(`/gene/${symbol}/`);
    },
    [router, setCommandOpen]
  );

  return (
    <Dialog
      open={commandOpen}
      onOpenChange={(open) => {
        setCommandOpen(open);
        if (!open) setQuery("");
      }}
    >
      <DialogHeader className="sr-only">
        <DialogTitle>Gene Search</DialogTitle>
        <DialogDescription>Search for genes by symbol or Ensembl ID</DialogDescription>
      </DialogHeader>
      <DialogContent
        className="top-1/3 translate-y-0 overflow-hidden rounded-xl! p-0 sm:max-w-lg"
        showCloseButton={false}
      >
        <Command shouldFilter={false}>
          <CommandInput
            placeholder="Search genes (symbol or Ensembl ID)..."
            value={query}
            onValueChange={setQuery}
          />
          <CommandList>
            <CommandEmpty>
              {query.length < 2
                ? "Type at least 2 characters..."
                : "No genes found."}
            </CommandEmpty>
            {results.length > 0 && (
              <CommandGroup heading="Genes">
                {results.map((gene) => (
                  <CommandItem
                    key={gene.symbol}
                    value={gene.symbol}
                    onSelect={() => handleSelect(gene.symbol)}
                    className="flex items-center justify-between"
                  >
                    <div className="flex items-center gap-2">
                      <span className="font-mono font-semibold">
                        {gene.symbol}
                      </span>
                      <span className="text-xs text-muted-foreground">
                        {gene.ensembl_id}
                      </span>
                    </div>
                    <div className="flex items-center gap-1">
                      {gene.is_deg && (
                        <Badge variant="secondary" className="text-xs">
                          DEG
                        </Badge>
                      )}
                      {gene.is_conserved_core && (
                        <Badge variant="secondary" className="text-xs">
                          CC
                        </Badge>
                      )}
                      <span className="text-xs text-muted-foreground">
                        {gene.layers_active || 0}/7
                      </span>
                    </div>
                  </CommandItem>
                ))}
              </CommandGroup>
            )}
            <CommandGroup heading="Navigation">
              <CommandItem
                onSelect={() => {
                  setCommandOpen(false);
                  router.push("/");
                }}
              >
                Home
              </CommandItem>
              <CommandItem
                onSelect={() => {
                  setCommandOpen(false);
                  router.push("/explore/");
                }}
              >
                Gene Explorer
              </CommandItem>
              <CommandItem
                onSelect={() => {
                  setCommandOpen(false);
                  router.push("/downloads/");
                }}
              >
                Downloads
              </CommandItem>
            </CommandGroup>
          </CommandList>
        </Command>
      </DialogContent>
    </Dialog>
  );
}
