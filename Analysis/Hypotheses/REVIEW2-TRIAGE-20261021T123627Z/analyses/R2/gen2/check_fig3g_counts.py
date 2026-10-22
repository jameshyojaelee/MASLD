# Blind check of gen2: split AN1 Fig 3G / universe transitions by the seven rerun-untestable genes.
# Reads only AN1 fig3g_states.tsv (3,776 genes with a finite best tier-1/2 SuSiE PP.H4 before the check).
import pandas as pd
D = "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/gene_level_20261021T165252Z/"
f = pd.read_csv(D + "fig3g_states.tsv", sep="\t", float_precision="round_trip")
print("genes", len(f))
print("rule a crosstab base x sens\n", pd.crosstab(f.state_a_base, f.state_a_sens))
print("rule b crosstab base x sens\n", pd.crosstab(f.state_b_base, f.state_b_sens))
x = f[f.finite_sens & (f.state_a_sens == "untestable")]
print("rule a: drawn and untestable after check", len(x), x.state.value_counts().to_dict(),
      "already untestable before", int((x.state_a_base == "untestable").sum()))
d = f[~f.finite_sens]
print("no finite PP.H4 after check", len(d), d.state.value_counts().to_dict(),
      "of which rerun-untestable (open circles)", int((d.state_a_base == "untestable").sum()))
