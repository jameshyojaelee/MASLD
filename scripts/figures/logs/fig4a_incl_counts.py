import os, sys
os.environ['FIG4A_INCLUSIVE']='1'
sys.path.insert(0, '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures')
import fig4a_v2_common as C
d, stats = C.get_data()
e = stats['enrich']
print("\n=== ENRICH FOLDS (inclusive) ===")
for grp in ('universe','coloc_only'):
    for lens in ('proteomics','spatial','scATAC'):
        x = e[grp][lens]
        print(f"  {grp:10s} {lens:10s} fold={x['fold']:.2f} p={x['p']:.1e} obs={x['obs']} exp={x['exp']:.1f}")
