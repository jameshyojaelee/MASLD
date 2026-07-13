import { create } from "zustand";
import type { AppState } from "./types";

export const useAppStore = create<AppState>((set) => ({
  // Transient default; the pre-hydration script + ThemeToggle set the real
  // theme (dark-first) from localStorage on mount.
  theme: "light",
  setTheme: (theme) => {
    if (typeof document !== "undefined") {
      document.documentElement.classList.toggle("dark", theme === "dark");
    }
    set({ theme });
  },

  sexLens: "combined",
  setSexLens: (sexLens) => set({ sexLens }),

  compareGenes: [],
  addCompareGene: (symbol) =>
    set((state) => {
      if (state.compareGenes.includes(symbol)) return state;
      if (state.compareGenes.length >= 5) return state;
      return { compareGenes: [...state.compareGenes, symbol] };
    }),
  removeCompareGene: (symbol) =>
    set((state) => ({
      compareGenes: state.compareGenes.filter((g) => g !== symbol),
    })),
  clearCompareGenes: () => set({ compareGenes: [] }),

  commandOpen: false,
  setCommandOpen: (commandOpen) => set({ commandOpen }),
}));
