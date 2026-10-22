# Scientific Figure Guidelines

Rules for generating publication-quality figures in this project. Style guidance is based on Fujii 2026 (Nature Human Behaviour, DOI: 10.1038/s41562-026-02466-9); output mechanics (paths, directory tree, naming) are project-specific. **All agents and scripts producing figures in this project MUST follow this document.**

> **Standalone Resource figure rules (2026-08-10):** Scientific identity and
> panel roles are governed by `docs/PAPER.md`. The five main
> figures are: (1) resource design/observability, (2) regulatory genetics, (3)
> cross-sectional established-state transcriptomics, (4) assay-native molecular
> and physical context, and (5) the MASLD Gene Catalog. Legacy
> directory and constant names such as `fig5_convergence` are compatibility
> paths only. They do not authorize a convergence ranking, causal-gene map,
> longitudinal progression, independent-validation funnel, discrete NMF
> subtype, or panel from the separate follow-up study. Existing PDFs are candidates or historical
> output files until a Plan 60 manifest accepts them.

For per-figure inventory and archive notes, see [`../figures/README.md`](../figures/README.md). This file is the single source of truth for figure rules; `figures/README.md` is the figure inventory.

## Core Principle

Every figure starts with a **key message**, not a statistical output. Before writing any plotting code, state in a comment: `# KEY MESSAGE: [one sentence describing what the reader should take away]`. If you cannot articulate the message, the figure is not ready to be made.

## Typography (MANDATORY — PI rule, enforced 2026-07-01)

**Every text element in every figure is font size 6, weight NORMAL — no bold, and no italic emphasis.** This is a hard rule and applies to ALL figures in this project: titles, axis labels, tick labels, legends, annotations, on-mark counts — everything is one size (6 pt). Do NOT use larger titles, bold headers/panel labels, or italics for emphasis. Set a single `FS = 6` constant (or `font.size: 6` in rcParams) and reuse it for every `text()`/label call; never pass `fontweight="bold"`. (Multi-panel letters A/B/C are added later in Illustrator during assembly, not baked into the script.)
> **Allowed exception — gene symbols may be italicized** (HGNC convention; `fontstyle="italic"` on the gene-name text only). Italics are for gene names, never for emphasis. Everything else stays roman. High-level overview panels carry NO gene names, so no italics appear there.

**Figures must be sized PROPORTIONAL to the font.** At font 6, a 7–8 inch-wide panel makes the text look tiny relative to the artwork — and it shrinks further during manuscript assembly, becoming illegible. Keep single panels compact, typically **~3–5.5 inches wide (≈75–140 mm)**, so 6 pt text reads at a sensible size relative to the whole panel. If a panel looks like "small font on a big canvas," **shrink the canvas — do not enlarge the font.**

Rationale: a manuscript figure is a composite of many panels; a single 6 pt standard + proportional canvas sizing keeps text uniform and legible across the whole figure instead of every panel having a different apparent text size. Reference implementation: `scripts/figures/fig4a_overview_candidates.py` (one `FS = 6`, no `fontweight`, compact `figsize`).

## Graph Type Selection

### Default replacements

| Instead of | Use | Why |
|------------|-----|-----|
| Bar chart of means | Dot/jitter/beeswarm plot with median | Bars hide distribution shape, sample size, outliers, modality |
| Bar chart of counts | Lollipop chart or Cleveland dot plot | Lower ink-to-data ratio, easier to compare across categories |
| Stacked bar chart (>4 groups) | Grouped dot plot or heatmap | Stacked bars beyond 3-4 groups are unreadable |
| Dynamite plot (bar + error bar) | Box + jitter, violin + points, or raincloud | Dynamite plots conceal everything about the distribution |
| Pie chart | Horizontal bar chart sorted by value | Humans are poor at judging angles; length is easier to compare |
| Dual y-axis line plot | Two aligned panels sharing the x-axis | Dual axes mislead by implying a relationship via scale choice |

### Bioinformatics-specific guidance

| Plot type | When to use | Key rules |
|-----------|-------------|-----------|
| Volcano plot | DEG results | Label only the genes that matter for your message. Do not label all significant genes. Use alpha (transparency) for non-significant points. Annotate thresholds directly on the plot, not only in the caption. |
| Heatmap | Expression matrices, correlation | Limit to top N genes/features relevant to the message. Always cluster rows AND columns unless order is biologically meaningful. Include a color scale bar with units. Never use rainbow colormap. |
| UMAP/tSNE | Dimensionality reduction | State perplexity/n_neighbors in the plot or caption. Use categorical colors for discrete groups. Do not overplot: subsample or use density contours if >50K points. Add a scale note that distances are not directly interpretable. |
| Dot plot (single-cell) | Marker gene expression | Limit to genes that support the message. Order cell types by biological hierarchy (not alphabetically). Size = percent expressing, color = mean expression. |
| Survival/KM curve | Time-to-event | Include number-at-risk table below the x-axis. Report log-rank p-value on the figure. Use distinct line styles in addition to color for accessibility. |
| Forest plot | Meta-analysis, effect sizes | Points + CI, never bars. Sort by effect size or by study date, not alphabetically. Include a vertical reference line at null effect. |
| Genome browser tracks | Epigenomic/regulatory data | Align all tracks to the same genomic coordinates. Label gene models clearly. Use consistent y-axis scales across comparable tracks. |
| Venn/Euler/UpSet | Set overlaps | Venn diagrams are acceptable for 2-3 sets only. For 4+ sets, use UpSet plots. Always report exact counts in each region, not just percentages. |

## Color Rules

### Choosing a scale type

| Data type | Scale | Example palettes |
|-----------|-------|-----------------|
| Categorical (cell types, conditions, clusters) | Qualitative | `Set2`, `Dark2`, `Paired`, `tab10` (limit to 8-10 colors max) |
| Sequential magnitude (expression level, p-value) | Sequential | `viridis`, `mako`, `cividis` (perceptually uniform, colorblind-safe) |
| Diverging from reference (log2FC, z-score) | Diverging | `RdBu_r`, `coolwarm`, `PiYG` (center at zero or reference value) |

### Mandatory constraints

- Never use rainbow/jet colormap. It is perceptually non-uniform and inaccessible.
- All figures must be interpretable under the three most common color vision deficiencies (protanopia, deuteranopia, tritanopia). Test with `colorspacious`, `colorblindr` (R), or https://www.color-blindness.com/coblis-color-blindness-simulator/.
- When using color to distinguish groups, also vary shape or linestyle as a redundant channel.
- Limit categorical colors to 8 per panel. If you have more groups, use faceting or a different encoding.
- For heatmaps, clip color scale to the 2nd-98th percentile range to prevent outliers from washing out the signal.

## Cognitive Load Reduction

### Remove

- Background shading and grid lines (use minimal or no grid)
- Redundant axis labels (if the title says "Gene Expression (TPM)", the y-axis label can be just "TPM")
- Box borders around legends
- 3D effects on any chart type
- Decorative elements that do not encode data

### Add

- Direct labels on lines/points instead of a distant legend (especially for <5 groups)
- Annotations pointing to the key comparison or result
- Panel labels (A, B, C) are added during Illustrator assembly, not in the script (see Typography rule — no bold text baked into panels)
- White space between panels (do not crowd subplots)

### Simplify

- One message per panel. If a panel makes two points, split it.
- If a figure has >6 panels, consider whether all are needed or if some belong in supplementary.
- Sort categories by the variable of interest (effect size, expression level), not alphabetically.

## Context-Dependent Complexity

Adjust figure detail based on where it will be shown:

| Context | Time | Interaction | Figure complexity |
|---------|------|-------------|-------------------|
| Lab meeting / group meeting | High | High | Raw data OK. Multiple panels OK. Can include QC plots and intermediate results. |
| Conference poster | Medium | Medium | Simplified. 1-2 key panels per figure. Verbal explanation fills gaps. |
| Oral presentation slide | Low | Low | Single message per slide. Large fonts. Annotate the takeaway directly on the figure. No dense heatmaps. |
| Journal article (main figure) | Self-paced | None | Moderate complexity. Clear structure. Every panel must be interpretable from the figure + caption alone. |
| Journal article (supplementary) | Self-paced | None | Higher complexity acceptable. QC, robustness checks, full parameter sweeps. |
| Twitter/social media | Very low | None | One panel. One message. Large text. Crop aggressively. |

## Tables vs Figures Decision

Use a **table** when:
- Reporting exact numerical values (model coefficients, sample sizes, p-values)
- Listing baseline characteristics or metadata
- Comparing <5 items across <5 dimensions
- The reader needs to look up specific values

Use a **figure** when:
- Showing distributions, trends, correlations, or spatial patterns
- Comparing across many groups or conditions
- The pattern is the message (not the exact numbers)
- Relationships between variables matter more than individual values

When unsure, ask: "Is the reader looking for a pattern or a number?" Pattern = figure. Number = table.

## Project-Specific Output Rules

These project rules override any generic guidance above when they conflict.

1. **PDF only. No PNG, no SVG.** Journal submission requires vector format; PDF is canonical here. (Supersedes the generic "save in PDF + PNG" advice.)
2. **All output paths MUST use the constants in [`scripts/figures/load_figure_data.R`](../scripts/figures/load_figure_data.R).** Source it at the top of every figure script. **Never hardcode a path.**
3. **Never write outputs to `scripts/figures/`.** Scripts live there; outputs go under `figures/`.
4. **Panels go flat in the parent figure directory.** No `panels/` subdirectory for final output.
5. **Do not create a new figure directory ad hoc.** Adding one requires (a) adding the constant + `dir.create` entry to `load_figure_data.R` and (b) updating the directory tree and path-constants table below.
6. **Do not modify figure scripts without explicit user request.** (Project standing rule.)

### Directory Tree

```
figures/
├── main/                              # Five main manuscript figures (reorg 2026-04-15: was six)
│   ├── fig1_atlas_overview/           # Fig 1 — Resource design, provenance, biological units, observability
│   ├── fig2_genetics/                 # Fig 2 — Phenotype-/ancestry-aware regulatory genetics
│   ├── fig3_RNAseq/                    # Fig 3 — Cross-sectional established-state remodeling + 117-program registry
│   ├── fig4_validation/               # Fig 4 — Assay-native protein, chromatin, and spatial/physical context
│   └── fig5_convergence/              # Fig 5 — MASLD Gene Catalog (legacy path name)
│
├── supplementary/                     # Supplementary figures
│   ├── figS01_qc_validation/          # S1  — QC, positive controls, variance partition
│   ├── figS02_progression/            # S2  — Cross-sectional stage contrasts and composition (legacy path name)
│   ├── figS03_deconvolution/          # S3  — Deconfounding (volcano C2 vs C13), deconvolution cross-method
│   ├── figS04_coloc/                  # S4  — COLOC extended, SuSiE comparison, cross-eQTL
│   ├── figS05_epigenomic_spatial/     # S5  — SCENIC+/chromVAR, scATAC, spatial transcriptomics
│   ├── figS06_cross_species/          # S6  — Historical/downstream cross-species context; not a human-paper claim
│   ├── figS07_ncrna/                  # S7  — ncRNA landscape, ceRNA, synteny
│   ├── figS08_subtyping_convergence/  # S8  — Continuous NMF axes + historical convergence sensitivities
│   ├── figS09_multi_ancestry/         # S9  — Multi-ancestry COLOC (FinnGen, BBJ)
│   ├── figS10_prediction/             # S10 — Historical predictive sensitivities; not a clinical-performance claim
│   ├── figS_sensitivity/              # Sens — DEG threshold sensitivity, ComBat-seq, age, BMI
│   ├── figS_hepatocyte_subtypes/      # Hepatocyte subclustering (UMAP, proportions, markers, spatial)
│   ├── figS_celltype_biology/         # Cell-type-resolved MASLD biology
│   ├── figS_therapeutics/             # Demoted Fig 6 panels — drug validation, THRB track, stratification
│   ├── figS_network/                  # v4 Architecture C: A–G architecture/replication/benchmark/hubs/neighborhoods/portal
│   ├── figS_rora_case_study/          # RORA multi-modal case study (Panels B–F; Panel A = fig3b locus zoom)
│   ├── figS_convergence_evidence/     # Historical heuristic-score diagnostics; never a main ranking
│   ├── figS_granular_staging/         # Cross-sectional granular-stage sensitivities (legacy filenames may say transition)
│   ├── figS_quantification/           # STAR vs Kallisto quantification comparison
│   └── archive/                       # Deprecated supplementary PDFs (provenance only)
│
├── misc/                              # Exploratory figures not intended for publication
├── archive/                           # Deprecated/provenance figures (do not delete)
│   ├── c2_deconv_raw_reference/       # Noncanonical raw-DEG C2 deconvolution diagnostics
│   └── rectangle_deconv_sensitivity/  # Matching Rectangle sensitivity diagnostics
├── examples/                          # Example panels for style reference
└── logs/                              # Figure generation logs
```

For per-figure inventory (panels, scripts, sources, archive subdirs) see [`../figures/README.md`](../figures/README.md).

### Path Constants Table

All constants are defined in `scripts/figures/load_figure_data.R` and resolve relative to `$MASLD_PROJECT_ROOT` (default `/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design`).

**CRITICAL: The canonical figure output root is `figures/` (NOT `figures_STAR/`).** `FIG_OUT` in `load_figure_data.R` is set to `figures/`. The `figures_STAR/` directory is a legacy leftover from the STAR-era pipeline and should NOT be used for new figures. All figure scripts MUST source `load_figure_data.R` and use its constants — never hardcode paths to either `figures/` or `figures_STAR/`.

| R Constant | Path (relative to project root) |
|------------|---------------------------------|
| `FIG_OUT` | `figures/` |
| `FIG_MAIN` | `figures/main/` |
| `FIG_SUPP` | `figures/supplementary/` |
| `FIG_MISC` | `figures/misc/` |
| `FIG_ARCHIVE` | `figures/archive/` |
| `FIG_DECONV_C2_ARCHIVE` | `figures/archive/c2_deconv_raw_reference/` |
| `FIG_DECONV_RECT_ARCHIVE` | `figures/archive/rectangle_deconv_sensitivity/` |
| `FIG1_DIR` | `figures/main/fig1_atlas_overview/` |
| `FIG2_DIR` | `figures/main/fig3_RNAseq/` (holds main **Fig 3** — constant name kept for back-compat) |
| `FIG3_BULK_DIR` | `figures/main/fig3_bulk_transcriptomics/` (current main **Fig 3** authority) |
| `FIG3_DIR` | `figures/main/fig2_genetics/` (holds main **Fig 2** — constant name kept for back-compat) |
| `FIG4_DIR` | `figures/main/fig4_validation/` |
| `FIG5_DIR` | `figures/main/fig5_convergence/` (holds main **Fig 5 MASLD Gene Catalog**; directory name is backward-compatible only) |
| `FIGS01_DIR` | `figures/supplementary/figS01_qc_validation/` |
| `FIGS02_DIR` | `figures/supplementary/figS02_progression/` (legacy path; cross-sectional stage analyses only) |
| `FIGS03_DIR` | `figures/supplementary/figS03_deconvolution/` |
| `FIGS04_DIR` | `figures/supplementary/figS04_coloc/` |
| `FIGS05_DIR` | `figures/supplementary/figS05_epigenomic_spatial/` |
| `FIGS06_DIR` | `figures/supplementary/figS06_cross_species/` |
| `FIGS07_DIR` | `figures/supplementary/figS07_ncrna/` |
| `FIGS08_DIR` | `figures/supplementary/figS08_subtyping_convergence/` |
| `FIGS09_DIR` | `figures/supplementary/figS09_multi_ancestry/` |
| `FIGS10_DIR` | `figures/supplementary/figS10_prediction/` (historical predictive sensitivities only) |
| `FIGS_SENS_DIR` | `figures/supplementary/figS_sensitivity/` |
| `FIGS_HEPSUB_DIR` | `figures/supplementary/figS_hepatocyte_subtypes/` |
| `FIGS_CELLTYPE_DIR` | `figures/supplementary/figS_celltype_biology/` |
| `FIGS_THERA_DIR` | `figures/supplementary/figS_therapeutics/` |
| `FIGS_NET_DIR` | `figures/supplementary/figS_network/` |
| `FIGS_RORA_DIR` | `figures/supplementary/figS_rora_case_study/` |
| `FIGS_CONV_EVID_DIR` | `figures/supplementary/figS_convergence_evidence/` |
| `FIGS_CAS13LIB_DIR` | `Cas13_Library_Design/figures/` (**separate future paper; excluded from the Resource figure manifest**) |
| `FIGS_GRANULAR_DIR` | `figures/supplementary/figS_granular_staging/` |
| `FIGS_QUANT_DIR` | `figures/supplementary/figS_quantification/` |

### Naming Conventions

| Figure type | Pattern | Example |
|-------------|---------|---------|
| Main figure panel | `fig{N}_{panel_letter}_{description}.pdf` | `fig1_a_cohort_overview.pdf` |
| Main figure composite | `fig{N}_composite.pdf` | `fig1_composite.pdf` |
| Supplementary panel | `figS{NN}_{description}.pdf` | `figS01_positive_controls.pdf` |
| Supplementary composite | `figS{NN}_composite.pdf` | `figS01_composite.pdf` |
| Sensitivity / misc | descriptive snake_case `.pdf` | `figS_deg_threshold_landscape.pdf` |

### How to Add a New Figure Directory

1. Add the constant to `scripts/figures/load_figure_data.R` and include it in the `dir.create` loop there.
2. Add a row to the path-constants table above and a line to the directory tree.
3. Source `scripts/figures/load_figure_data.R` in your figure script and use the new constant — never the literal path.
4. Write output: `ggsave(file.path(NEW_DIR, "figSXX_description.pdf"), ...)`.

## Implementation Checklist

Before finalizing any figure, verify:

- [ ] Panel appears in the Plan 60 accepted-panel manifest for the correct Resource role
- [ ] Cross-sectional stage data are not labeled longitudinal progression, transition, cascade, or treatment window
- [ ] COLOC is labeled regulatory evidence/candidate nomination, not causal-gene proof
- [ ] Assay-native context is not called independent or functional validation without source- and sample-level justification
- [ ] Figure 5 defaults to evidence state, testability, source dependence, and next experiment—not score or rank
- [ ] No separate follow-up study roster, assay design, phenotype, hit, or validation claim appears in the Resource figure set
- [ ] Key message is stated and the figure supports it
- [ ] Graph type matches the data structure (not just habit)
- [ ] Color scheme matches data type (categorical/sequential/diverging)
- [ ] Colorblind accessible (tested or using known-safe palette)
- [ ] No rainbow/jet colormap
- [ ] No bar charts hiding distributions (unless counts/proportions)
- [ ] Axis labels have units
- [ ] ALL text is font size 6, weight normal (no bold, no italics); canvas is proportional to the font (compact ~3–5.5 in wide, not tiny text on a big canvas)
- [ ] Legend is close to the data it describes (or replaced with direct labels)
- [ ] Non-essential gridlines, borders, and decorations removed
- [ ] Panels are labeled (A, B, C) if multi-panel
- [ ] Figure resolution is >= 300 DPI for print, >=150 DPI for screen
- [ ] Saved as **PDF only** (vector); no PNG/SVG companions
- [ ] Output path uses a constant from `load_figure_data.R` (no hardcoded paths)
- [ ] Output is under `figures/`, not `scripts/figures/`

## Python Defaults

When generating figures with matplotlib/seaborn, apply these defaults at the top of every plotting script:

```python
import matplotlib.pyplot as plt
import seaborn as sns

# -- Figure defaults --
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.grid': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
    # Typography rule: ONE size (6 pt), no bold/italic anywhere. See "Typography" above.
    'font.size': 6,
    'axes.titlesize': 6,
    'axes.labelsize': 6,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'axes.titleweight': 'normal',
    'axes.labelweight': 'normal',
    'pdf.fonttype': 42,      # editable text in Illustrator
    'ps.fonttype': 42,
    'legend.frameon': False,
    # Keep the canvas PROPORTIONAL to the 6 pt font — compact single panels
    # (~3–5.5 in wide). Do not leave a big canvas with tiny text.
    'figure.figsize': (4.5, 3.0),
})

# -- Default palette (colorblind-safe) --
sns.set_palette('Set2')
```

## R Defaults

When generating figures with ggplot2, apply this theme:

```r
library(ggplot2)

# Typography rule: base_size 6, weight NORMAL (no bold/italic), canvas proportional.
theme_publication <- theme_minimal(base_size = 6) +
  theme(
    text = element_text(size = 6, face = "plain", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 6),
    plot.title = element_text(face = "plain", size = 6),   # NOT bold
    plot.margin = margin(6, 6, 6, 6)
  )
# Save compact: ggsave(..., width = 4.5, height = 3.0, units = "in")  # canvas proportional to 6 pt font

theme_set(theme_publication)

# Colorblind-safe discrete palette
scale_colour_cb <- scale_colour_brewer(palette = "Set2")
scale_fill_cb <- scale_fill_brewer(palette = "Set2")
```
