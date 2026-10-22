# Pilot: clear-vacuole (macrovesicular fat) area fraction on GTEx liver H&E, ~2 um/px level.
# Parameters fixed before any result was seen:
#   tissue: HSV saturation>=20/255 at 4x-downsampled grid, holes filled, objects < 2000 px (4x grid) removed
#   white: V>=200 and S<=40 (0-255) inside filled tissue mask
#   droplet: 8-connected white component, area 20..1300 px at ~1.98 um/px (diameter ~10-80 um),
#            bbox extent >= 0.55, bbox aspect <= 2.5
#   fat fraction = droplet area / tissue area
import sys, time, os, numpy as np, pandas as pd, cv2, wsidicom
from scipy import ndimage as ndi
from multiprocessing import Pool
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT=os.path.join(os.getcwd(),'morph_fat.tsv')
T0=time.time(); BUDGET=float(sys.argv[1]) if len(sys.argv)>1 else 380
def one(args):
    donor, uid = args
    if time.time()-T0 > BUDGET: return dict(donor=donor,uid=uid,status='skipped_time')
    try:
        w=wsidicom.WsiDicom.open(B+'data/external/gtex_liver_slides/series/'+uid)
        lv=min(w.levels, key=lambda l: abs(l.mpp.width-2.0))
        a=np.asarray(w.read_region((0,0),lv.level,(lv.size.width,lv.size.height)).convert('RGB'))
        mx=a.max(2).astype(np.int16); mn=a.min(2).astype(np.int16)
        S=np.where(mx>0,(mx-mn)*255//np.maximum(mx,1),0).astype(np.uint8)
        small=cv2.resize(S,(S.shape[1]//4,S.shape[0]//4),interpolation=cv2.INTER_AREA)
        tis=small>=20
        tis=ndi.binary_fill_holes(tis)
        lab,nl=ndi.label(tis); sizes=np.bincount(lab.ravel()); keep=sizes>=2000; keep[0]=False
        tis=keep[lab]
        tisF=cv2.resize(tis.astype(np.uint8),(S.shape[1],S.shape[0]),interpolation=cv2.INTER_NEAREST).astype(bool)
        white=(mx>=200)&(S<=40)&tisF
        n,lab2,st,_=cv2.connectedComponentsWithStats(white.astype(np.uint8),connectivity=8)
        st=st[1:]; area=st[:,4]; wd=st[:,2]; ht=st[:,3]
        ext=area/(wd*ht); asp=np.maximum(wd,ht)/np.maximum(np.minimum(wd,ht),1)
        dro=(area>=20)&(area<=1300)&(ext>=0.55)&(asp<=2.5)
        tarea=int(tisF.sum())
        return dict(donor=donor,uid=uid,status='ok',mpp=round(lv.mpp.width,3),tissue_px=tarea,
                    white_frac=float(white.sum()/max(tarea,1)),fat_frac=float(area[dro].sum()/max(tarea,1)),
                    n_droplets=int(dro.sum()),secs=round(time.time()-T0,1))
    except Exception as e:
        return dict(donor=donor,uid=uid,status='error:'+repr(e)[:200])
if __name__=='__main__':
    order=pd.read_csv('order.tsv',sep='\t')
    rows=[]
    with Pool(4) as p:
        for r in p.imap(one, list(zip(order.donor,order.uid))):
            rows.append(r)
            if len(rows)%10==0: pd.DataFrame(rows).to_csv(OUT,sep='\t',index=False)
    pd.DataFrame(rows).to_csv(OUT,sep='\t',index=False)
    print('done',len(rows),'ok',sum(r['status']=='ok' for r in rows),'elapsed',round(time.time()-T0,1))
