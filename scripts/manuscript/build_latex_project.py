#!/usr/bin/env python3
"""Prepare the editable manuscript and current PDF panels for Overleaf."""
from pathlib import Path
import csv
import hashlib
import json
import re
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'docs/manuscript/latex'
PANDOC = Path('/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pandoc')
DRAFT = ROOT / 'docs/manuscript/draft'


def convert(text):
    result = subprocess.run([str(PANDOC), '-f', 'markdown+autolink_bare_uris', '-t', 'latex',
                           '--wrap=preserve'], input=text, text=True,
                          capture_output=True, check=True).stdout
    result = result.replace('\\_', '\\_\\allowbreak{}')
    result = result.replace('AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3',
        'AllPredictions.\\allowbreak{}AvgHiC.\\allowbreak{}ABC0.015.\\allowbreak{}minus150.\\allowbreak{}ForABCPaperV3')
    return result


def write(path, text):
    dest = OUT / path
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text)


def main():
    if OUT.exists():
        raise SystemExit(f'Refusing to overwrite existing project: {OUT}')
    OUT.mkdir(parents=True)
    sources = ['abstract', '01_intro', 'fig1', 'fig2', 'fig3', 'fig4', 'fig5',
               'fig5_discussion', 'fig7', 'references', 'supplementary_notes']
    for name in sources:
        write(f'source/{name}.md', (DRAFT / f'{name}.md').read_text())
    for name in ['METHODS', '05_figure_legends']:
        write(f'source/{name}.md', (ROOT / f'docs/manuscript/{name}.md').read_text())
    catalog, discussion = (DRAFT / 'fig5_discussion.md').read_text().split('## Discussion\n', 1)
    texts = {name: (DRAFT / f'{name}.md').read_text() for name in sources}
    texts['catalog'] = catalog
    texts['discussion'] = '# Discussion\n\n' + discussion
    texts['methods'] = (ROOT / 'docs/manuscript/METHODS.md').read_text()
    texts['legends'] = (ROOT / 'docs/manuscript/05_figure_legends.md').read_text()
    order = ['abstract', '01_intro', 'fig1', 'fig2', 'fig3', 'fig4', 'fig5',
             'catalog', 'fig7', 'discussion', 'references', 'methods',
             'supplementary_notes', 'legends']
    for name in order:
        write(f'sections/{name}.tex', convert(texts[name]))
    # Preserve the existing curated bibliography verbatim in the compiled text.
    # BibTeX records provide editable DOI/note fields without inventing metadata.
    refs = texts['references'].strip().split('\n\n')[1:]
    bib = []
    for n, ref in enumerate(refs, 1):
        doi = re.search(r'https://doi.org/(\S+)', ref)
        note = re.sub(r'^\d+\.\s*', '', ref).replace('*', '')
        note = convert(note).strip()
        fields = [f'  note = {{{note}}}']
        if doi:
            fields.append(f'  doi = {{{doi[1].rstrip(".")}}}')
        bib.append('@misc{resource_ref_%02d,\n%s\n}' % (n, ',\n'.join(fields)))
    write('references.bib', '\n\n'.join(bib) + '\n')

    main_panels = {}
    f = ROOT / 'figures/main'
    main_panels[1] = [
        ('A', ROOT/'figures/candidates/fig1a-resource-proposal-22238205/fig1a_resource_schematic.pdf', 'proposal'),
        ('B', f/'fig1_atlas_overview/fig1b_complementary_maps.pdf', 'working'),
        ('C', f/'fig1_atlas_overview/fig1c_genetic_trait_scope.pdf', 'working'),
        ('D', ROOT/'figures/candidates/fig1-review-proposal-22238377/fig1d_analysis_workflow.pdf', 'proposal')]
    names2 = ['fig2A_gwas_cascade', 'Fig2B_PIP_vs_SuSiE-coloc',
              'Fig2C_coding_noncoding_split', 'Fig2D_pip_architecture_by_trait_directness',
              'Fig2E_ancestry_unique_coloc_GWS', 'Fig2F_crossancestry_coloc',
              'Fig2G_RORA_NAFLD', 'Fig2H_FABP1_ALT', 'Fig2I_locus_legend']
    main_panels[2] = [(chr(65+i), f/'fig2_genetics/panels'/f'{n}.pdf', 'working') for i,n in enumerate(names2)]
    main_panels[2].insert(4, ('D ancestry', f/'fig2_genetics/panels/Fig2D_coding_by_ancestry.pdf', 'working_companion'))
    names3 = ['fig3a_cohort_metadata_matrix', 'fig3b_nas_fib_grid', 'fig3c_cohort_alluvial',
              'fig3d_pca_fibrosis_gradient', 'fig3e_stage_remodeling',
              'fig3f_stage_deg_genetics_matrix', 'fig3g_genetic_state_complementarity']
    main_panels[3] = [(chr(65+i), f/'fig3_bulk_transcriptomics/panels'/f'{n}.pdf', 'working') for i,n in enumerate(names3)]
    main_panels[3].insert(5, ('E replication', f/'fig3_bulk_transcriptomics/panels/fig3e_steatosite_replication_strip.pdf', 'working_companion'))
    main_panels[4] = [(chr(65+i), p, 'historical_atlas_input' if i<4 else 'working') for i,p in enumerate(sorted((f/'fig4_singlecell_programs/panels').glob('fig4[a-h]_*.pdf')))]
    main_panels[5] = [(chr(65+i), p, 'working') for i,p in enumerate(sorted((f/'fig5_molecular_context/panels').glob('fig5[a-d]_*.pdf')))]
    main_panels[6] = [(chr(65+i), p, 'working') for i,p in enumerate(sorted((f/'fig6_gene_catalog/panels').glob('fig6[a-e]_*.pdf')))]
    main_panels[7] = [(chr(65+i), p, 'working') for i,p in enumerate(sorted((f/'fig7_model_framework').glob('fig7[A-F]_*.pdf')))]
    expected = {1:4,2:10,3:8,4:8,5:4,6:5,7:6}
    assert {k:len(v) for k,v in main_panels.items()} == expected
    manifest = []

    def copy_pdf(p, dest, role, state):
        if not p.is_file():
            raise FileNotFoundError(p)
        target = OUT / dest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, target)
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        manifest.append({'file':dest,'source':str(p.relative_to(ROOT)),
                         'role':role,'state':state,'sha256':digest})

    # LaTeX review layouts use current individual PDFs rather than stale composites.
    # Original panel PDFs remain separately editable in Illustrator.
    plates = ['\\clearpage\n\\section{Main figures}\n']
    for num, panels in main_panels.items():
        page_break = '' if num == 1 else '\\clearpage\n'
        plates += [page_break + f'\\subsection{{Figure {num}}}\n',
                   '\\noindent Working panel layout; full captions are in Figure legends.\\par\n']
        if num == 1:
            plates += ['\\noindent Panels A and D remain overview proposals.\\par\n']
        if num == 4:
            plates += ['\\noindent Atlas-derived panels retain the repeated-capture qualification.\\par\n']
        for label,p,state in panels:
            dest = f'figures/main/figure{num}/{p.name}'
            copy_pdf(p,dest,f'Figure {num}{label}',state)
            plates += ['\\noindent\\begin{minipage}{\\linewidth}\n',
                       f'\\textbf{{{label}}}\\par\n',
                       '\\begin{center}\n',
                       f'\\includegraphics[width=0.96\\linewidth,height=0.68\\textheight,keepaspectratio]{{{dest}}}\n',
                       '\\end{center}\n\\end{minipage}\\par\\medskip\n']
    write('sections/main_figures.tex', ''.join(plates))

    supplements = []
    supplements += [p for p in (f/'fig2_genetics/panels').glob('FigS2*.pdf') if not p.name.startswith(('FigS2E_', 'FigS2N_'))]
    supplements += list((f/'fig3_bulk_transcriptomics/panels/supplementary').glob('*.pdf'))
    supplements += list((f/'fig3_bulk_transcriptomics/panels').glob('figs3*.pdf'))
    supplements += list((f/'fig4_singlecell_programs/panels/supplementary').glob('*.pdf'))
    supplements += list((f/'fig5_molecular_context/panels').glob('figS5*.pdf'))
    supplements += list((f/'fig5_molecular_context/panels/_archive/demoted-from-main-20261103').glob('*.pdf'))
    suptex = ['\\clearpage\n\\section{Supplementary figure artwork}\n',
              'Panel filenames identify the working artwork. Unassigned letters and historical-input qualifications remain as stated in the legends.\\par\n']
    for p in sorted(set(supplements)):
        dest = 'figures/supplementary/' + str(p.relative_to(f))
        copy_pdf(p,dest,'Supplementary artwork','working_or_qualified_as_in_legend')
        name = p.stem.replace('_', ' ')
        suptex += ['\\clearpage\n', f'\\subsection*{{{name}}}\n',
                   '\\begin{center}\n',
                   f'\\includegraphics[width=\\linewidth,height=0.83\\textheight,keepaspectratio]{{{dest}}}\n',
                   '\\end{center}\n']
    write('sections/supplementary_figures.tex', ''.join(suptex))
    header = '''\n\\usepackage{graphicx}
\\setmainfont[Path=fonts/,BoldFont=DejaVuSerif-Bold.ttf,ItalicFont=DejaVuSerif.ttf,ItalicFeatures={FakeSlant=0.2},BoldItalicFont=DejaVuSerif-Bold.ttf,BoldItalicFeatures={FakeSlant=0.2}]{DejaVuSerif.ttf}
\\usepackage{longtable,booktabs,array,calc}
\\usepackage{etoolbox}
\\usepackage{fancyvrb}
\\newcommand{\\passthrough}[1]{#1}
\\newcommand{\\pandocbounded}[1]{#1}
\\AtBeginEnvironment{longtable}{\\footnotesize}
\\setlength{\\emergencystretch}{4em}
\\sloppy
'''
    write('preamble-extra.tex', header)
    for name in ['DejaVuSerif.ttf', 'DejaVuSerif-Bold.ttf']:
        dest = OUT / 'fonts' / name
        dest.parent.mkdir(exist_ok=True)
        shutil.copyfile(Path('/usr/share/fonts/truetype/dejavu')/name, dest)
    shutil.copyfile('/usr/share/doc/fonts-dejavu-core/copyright', OUT/'fonts/LICENSE.txt')
    title = '''# A multimodal human MASLD resource maps inherited genetic variation and multicellular disease remodeling

James Lee^1,2,\\*^, Victor Passanisi^1,2,\\*^, Neha Saravanan^1,2^, Sofia Gonzalez-Nolde^1,2^, Neville E. Sanjana^1,2,\\#^

^1^ New York Genome Center, New York, NY  
^2^ Department of Biology, New York University, New York, NY

\\* Equal contribution  
\\# Correspondence: nsanjana@nygenome.org (NES)

'''
    body = convert(title)
    body += '\n'.join(f'\\input{{sections/{n}.tex}}' for n in order)
    body += '\n\\input{sections/main_figures.tex}\n\\input{sections/supplementary_figures.tex}\n'
    template = subprocess.run([str(PANDOC), '-f','markdown','-t','latex','--standalone',
        '--pdf-engine=xelatex','-V','geometry:margin=1in',
        '--include-in-header',str(OUT/'preamble-extra.tex')], input='',text=True,
        capture_output=True,check=True).stdout
    template = template.replace('\\begin{document}', '\\begin{document}\n'+body)
    write('main.tex', template)
    write('latexmkrc', "$pdflatex = 'xelatex %O %S';\n$pdf_mode = 1;\n")
    with (OUT/'figure_inventory.tsv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
    write('README.txt', '''Main document: main.tex. Compiler: XeLaTeX.
Edit sections/*.tex directly. source/ contains the matching Markdown snapshot.
Main figure PDFs are in figures/main/figure1 through figure7. Original PDF panels
are preserved; LaTeX controls their placement. Supplementary artwork is included.
Figure 1 A/D remain proposals. Genetic-dependent claims remain provisional and
atlas-derived panels retain the repeated-physical-capture qualification.
References are rendered from the existing curated bibliography; references.bib
provides DOI/note records for citation editing. Citation placeholders, unassigned
panel letters, deposit accessions and author-completion fields remain unresolved.
This is an editable working manuscript, not a new adopted scientific release.
''')
    checks={'main_figures':{str(k):len(v) for k,v in main_panels.items()},
            'pdf_files':len(manifest),'supplementary_pdf_files':len(set(supplements)),
            'source_snapshot_files':len(sources)+2}
    write('project_checks.json',json.dumps(checks,indent=2)+'\n')
    archive=OUT.parent/'masld-overleaf-20261108.zip'
    if archive.exists():
        raise FileExistsError(archive)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in OUT.rglob('*'):
            if p.is_file():z.write(p,p.relative_to(OUT))
    print(json.dumps({'project':str(OUT),'archive':str(archive),**checks},indent=2))


if __name__ == '__main__':
    main()
