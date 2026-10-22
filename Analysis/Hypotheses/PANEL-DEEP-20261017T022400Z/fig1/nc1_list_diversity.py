"""How different are the expression-matched NC-1 style lists from each other (seeded)."""
import numpy as np, pandas as pd, json
ev = pd.read_csv("recompute_1b_U_eval_genes.tsv.gz", sep="\t")
gen = ev[ev.is_gen]; non = ev[~ev.is_gen].reset_index(drop=True)
rng = np.random.default_rng(20261017)
A = non.AveExpr.to_numpy(); lists = []
for i in range(200):
    used = np.zeros(len(non), bool); pick = []
    for a in gen.AveExpr.to_numpy()[rng.permutation(len(gen))]:
        d = np.abs(A - a); d[used] = np.inf; j = int(np.argmin(d)); used[j] = True; pick.append(j)
    lists.append(frozenset(pick))
union = frozenset().union(*lists); inter = frozenset.intersection(*lists)
degs_union = int(non.loc[list(union), "is_deg"].sum())
res = dict(n_lists=200, list_size=426, union_genes=len(union), intersection_genes=len(inter),
           distinct_lists=len(set(lists)), deg_in_intersection=int(non.loc[list(inter), "is_deg"].sum()), deg_in_union=degs_union,
           smd_aveexpr=float((non.loc[list(lists[0]), "AveExpr"].mean() - gen.AveExpr.mean()) / gen.AveExpr.std()))
json.dump(res, open("nc1_list_diversity.json", "w"), indent=1); print(res)
