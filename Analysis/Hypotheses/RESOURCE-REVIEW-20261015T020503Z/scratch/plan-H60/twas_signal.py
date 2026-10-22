# TWAS-side only. Does not touch DE t. Counts signal per trait and file mtimes.
import glob, os, numpy as np, pandas as pd, time
OT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/causal_inference/otters_broadaway/"
for t in ["ukbb_alt","ghodsian","ukbb_ast","ukbb_ggt","pdff"]:
    fs=sorted(glob.glob(OT+f"{t}/chr*/P0.05.txt"))
    if not fs: print(t,"no P0.05 files"); continue
    d=pd.concat([pd.read_csv(f,sep="\t") for f in fs]).dropna(subset=["FUSION_Z"])
    mt=max(os.path.getmtime(f) for f in fs)
    print(f"{t}: files {len(fs)} genes {d.TargetID.nunique()} mean z^2 {np.mean(d.FUSION_Z**2):.3f} |z|>3 {int((d.FUSION_Z.abs()>3).sum())} expected {len(d)*2*0.0013499:.0f} latest mtime {time.strftime('%Y-%m-%d', time.localtime(mt))}")
