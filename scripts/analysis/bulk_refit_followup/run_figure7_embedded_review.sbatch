#!/bin/bash
#SBATCH --job-name=figure7_embed_review
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/figure7_embed_review_%j.out
#SBATCH --error=logs/figure7_embed_review_%j.err
set -euo pipefail
export FIG34_OUT="$PWD/figures/candidates/figure7-embedded-native-review-proposal-${SLURM_JOB_ID}"
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
mkdir -p "$FIG34_OUT/source"
printf '%s\n' "$FIG34_OUT" > "logs/figure7_embed_review_${SLURM_JOB_ID}.root"
cp scripts/analysis/bulk_refit_followup/assemble_figure7_review.py scripts/analysis/bulk_refit_followup/run_figure7_embedded_review.sbatch "$FIG34_OUT/source/"
python "$FIG34_OUT/source/assemble_figure7_review.py"
printf '/Helvetica (/gpfs/commons/home/jameslee/.local/share/fonts/Helvetica.ttf) ;\n' > "$FIG34_OUT/Fontmap"
/usr/bin/gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -dEmbedAllFonts=true -dSubsetFonts=true -dPDFSETTINGS=/prepress -sFONTMAP="$FIG34_OUT/Fontmap" -sOutputFile="$FIG34_OUT/Figure7_embedded_native_review_PROPOSAL.pdf" "$FIG34_OUT/Figure7_review_proposal.pdf"
python - <<'PY'
import importlib.util,json,os,subprocess,hashlib,xml.etree.ElementTree as ET
from collections import Counter
import numpy as np
from scipy.optimize import linear_sum_assignment
from pathlib import Path
from pypdf import PdfReader
out=Path(os.environ['FIG34_OUT']);old=out/'Figure7_review_proposal.pdf';new=out/'Figure7_embedded_native_review_PROPOSAL.pdf'
exe='/gpfs/commons/home/jameslee/.local/bin/'
def boxes(path):
    root=ET.fromstring(subprocess.check_output([exe+'pdftotext','-bbox',str(path),'-']))
    return [(e.text,tuple(float(e.attrib[k]) for k in ['xMin','yMin','xMax','yMax'])) for e in root.iter() if e.tag.endswith('word')]
a,b=boxes(old),boxes(new);assert Counter(x[0] for x in a)==Counter(x[0] for x in b)
# PDF conversion can change extraction order of adjacent/rotated words.
# Match equal words by their page positions rather than extraction sequence.
shift=0.; added_letter_bbox_shift=0.; added_letter_count=0
for word in set(x[0] for x in a):
    aa=np.array([x[1] for x in a if x[0]==word]);bb=np.array([x[1] for x in b if x[0]==word])
    cost=np.max(np.abs(aa[:,None,:]-bb[None,:,:]),axis=2)
    ii,jj=linear_sum_assignment(cost)
    for i,j in zip(ii,jj):
        # Six added panel letters are outside the native panel rectangles.
        # GS changes their embedded-font descriptor bbox, not source artwork.
        is_added_letter=word in 'ABCDEF' and min(abs(aa[i,0]-x) for x in [18.,416.])<.01 and min(abs(aa[i,1]-y) for y in [5.670875,227.670875,477.670875])<.01
        if is_added_letter:
            added_letter_count+=1;added_letter_bbox_shift=max(added_letter_bbox_shift,float(cost[i,j]))
        else: shift=max(shift,float(cost[i,j]))
assert added_letter_count==6
assert shift<=.25,shift
assert list(PdfReader(old).pages[0].mediabox)==list(PdfReader(new).pages[0].mediabox)
spec=importlib.util.spec_from_file_location('assembler','figures/candidates/fig3-fig4-assembled-774-20261013T074505Z/producers/assemble_fig34.py');asm=importlib.util.module_from_spec(spec);spec.loader.exec_module(asm)
fonts=asm.fonts_of(PdfReader(new));assert all(fonts.values()) and all(k[1]!='/Type3' for k in fonts)
check={'text_identical':True,'page_box_identical':True,'native_panel_scale':1,'max_native_word_bbox_difference_pt':shift,'added_letter_bbox_descriptor_difference_pt':added_letter_bbox_shift,'added_letter_count':added_letter_count,'all_fonts_embedded':True,'no_type3_fonts':True,'fonts':[{'font':k[0],'subtype':k[1],'embedded':v} for k,v in fonts.items()],'sha256':hashlib.sha256(new.read_bytes()).hexdigest(),'proposal_only':True}
(out/'embedded_identity_font_checks.json').write_text(json.dumps(check,indent=2)+'\n')
subprocess.run([exe+'pdftoppm','-scale-to','1600','-singlefile','-png',str(new),str(out/'proofs/Figure7_embedded_native_review_PROPOSAL')],check=True)
print(json.dumps(check,indent=2))
PY
