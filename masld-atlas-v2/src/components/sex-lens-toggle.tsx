"use client";

import { useAppStore } from "@/lib/store";

const OPTIONS = [
  { value: "female" as const, label: "F", activeColor: "text-pink-600" },
  { value: "combined" as const, label: "All", activeColor: "" },
  { value: "male" as const, label: "M", activeColor: "text-blue-600" },
] as const;

export function SexLensToggle() {
  const { sexLens, setSexLens } = useAppStore();

  return (
    <div className="flex items-center gap-2">
      <span className="text-[10px] font-semibold uppercase tracking-widest text-muted-foreground">
        Sex
      </span>
      <div className="flex rounded-md border border-border overflow-hidden">
        {OPTIONS.map((opt) => {
          const isActive = sexLens === opt.value;
          return (
            <button
              key={opt.value}
              onClick={() => setSexLens(opt.value)}
              className={[
                "px-2 py-1 text-xs font-medium transition-colors",
                isActive
                  ? `bg-primary text-primary-foreground ${opt.activeColor}`
                  : "text-muted-foreground hover:text-foreground hover:bg-muted",
              ].join(" ")}
            >
              {opt.label}
            </button>
          );
        })}
      </div>
    </div>
  );
}
