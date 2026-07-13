"use client";

import { useEffect, useRef } from "react";
import { Button } from "@/components/ui/button";
import { useAppStore } from "@/lib/store";

export function ThemeToggle() {
  const { theme, setTheme } = useAppStore();

  useEffect(() => {
    // Initialize the store from the theme the pre-hydration script already
    // applied (localStorage, else the class it set from the OS preference),
    // so a first-visit OS-dark preference is not reset to light.
    const saved = localStorage.getItem("masld-atlas-theme");
    const initial =
      saved === "dark" || saved === "light"
        ? saved
        : document.documentElement.classList.contains("dark")
          ? "dark"
          : "light";
    setTheme(initial);
  }, [setTheme]);

  // Persist + sync the DOM only on actual theme CHANGES, never on mount — the
  // pre-hydration script already applied the correct class, and running this on
  // mount (with the store's transient initial value) would clobber the user's
  // stored choice.
  const firstRun = useRef(true);
  useEffect(() => {
    if (firstRun.current) {
      firstRun.current = false;
      return;
    }
    localStorage.setItem("masld-atlas-theme", theme);
    document.documentElement.classList.toggle("dark", theme === "dark");
  }, [theme]);

  return (
    <Button
      variant="ghost"
      size="sm"
      onClick={() => setTheme(theme === "light" ? "dark" : "light")}
      className="h-8 w-8 px-0"
      aria-label={`Switch to ${theme === "light" ? "dark" : "light"} mode`}
    >
      {theme === "light" ? (
        <svg
          className="h-4 w-4"
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
          strokeWidth={2}
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            d="M20.354 15.354A9 9 0 018.646 3.646 9.003 9.003 0 0012 21a9.003 9.003 0 008.354-5.646z"
          />
        </svg>
      ) : (
        <svg
          className="h-4 w-4"
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
          strokeWidth={2}
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            d="M12 3v1m0 16v1m9-9h-1M4 12H3m15.364 6.364l-.707-.707M6.343 6.343l-.707-.707m12.728 0l-.707.707M6.343 17.657l-.707.707M16 12a4 4 0 11-8 0 4 4 0 018 0z"
          />
        </svg>
      )}
    </Button>
  );
}
