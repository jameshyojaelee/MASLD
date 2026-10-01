"""Export permitted development beta only after fixed metadata admission.

No mixed-fold source, p/q values, old predictions, fitting or evaluation.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT/'GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/allele_stack_inputs_21849187/development_labels_folds1to4.tsv.gz'
FIXED = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/target_region_manifest_21995444'
MANIFEST = FIXED/'development_target_manifest_1mb.tsv.gz'
EXCLUDED = FIXED/'geometry_exclusions.tsv'
HASHES = {'source':'72bfc08cb05dc0ff8d4e69580b1d20484ea6f9f2a3bb45336c8336edea1e9916',
          'manifest':'ef47cc5b60aa2aab542c686772a620575af3f90e1567eac7c17beb1ac2c63eef',
          'geometry_exclusions':'be2b83270a81650f0e8764dda93f85b688b5085c230018f5d06669bbc03f6f80'}
SOURCE_META = ['lead_variant_id','chr','pos_hg38','ref','alt','peak_id',
               'peak_start_hg38','peak_stop_hg38','heldout_fold','block_1mb']
GEOMETRY_META = ['variant_id','peak_id','chromosome','heldout_fold','length_bp',
                 'peak_start0','peak_end0','coordinate_state','target_readout_state']
OUTPUT = ['key','lead_variant_id','chr','pos_hg38','ref','alt','peak_id',
          'peak_start0','peak_end0','heldout_fold','block_1mb','beta_alt']


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def require(value, message):
    if not value:
        raise ValueError(message)


def integer(frame, names):
    for name in names:
        require(frame[name].str.fullmatch(r'[0-9]+').all(), 'Noninteger metadata: '+name)
        frame[name] = pd.to_numeric(frame[name], errors='raise').astype(np.int64)


def main():
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Compute allocation required.')
    paths = {'source':SOURCE, 'manifest':MANIFEST, 'geometry_exclusions':EXCLUDED}
    for name, path in paths.items():
        require(sha(path) == HASHES[name], 'Fixed input changed: '+name)
    # First pass imports identity metadata only. No beta or significance field.
    source = pd.read_csv(SOURCE, sep='\t', usecols=SOURCE_META, dtype=str, keep_default_na=False)
    integer(source, ['pos_hg38','peak_start_hg38','peak_stop_hg38','heldout_fold'])
    require(len(source) == 25748, 'Development source must have exactly 25748 rows.')
    require(set(source.heldout_fold) == {1,2,3,4}, 'Forbidden/missing development fold.')
    source['chr'] = 'chr'+source['chr'].str.removeprefix('chr')
    require(source['chr'].str.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])').all(), 'Nonautosomal chromosome.')
    require(source.pos_hg38.gt(0).all(), 'Variant position must be one-based positive.')
    require(source.block_1mb.str.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2]):[0-9]+').all(), 'Invalid source block.')
    require(source.block_1mb.str.split(':').str[0].str.removeprefix('chr').equals(
            source['chr'].str.removeprefix('chr')), 'Source block chromosome mismatch.')
    require(source.ref.str.fullmatch('[ACGT]').all() and source.alt.str.fullmatch('[ACGT]').all()
            and source.ref.ne(source.alt).all(), 'Invalid source SNV alleles.')
    source['key'] = (source['chr'].str.removeprefix('chr')+':'+source.pos_hg38.astype(str)
                     +':'+source.ref+':'+source.alt)
    require(source.lead_variant_id.str.removeprefix('chr').equals(source.key), 'Variant ID differs from metadata.')
    require(not source.key.duplicated().any() and not source.lead_variant_id.duplicated().any(), 'Repeated variant key.')
    require(source.groupby('chr').heldout_fold.nunique().max() == 1, 'Chromosome crosses folds.')
    source = source.rename(columns={'peak_start_hg38':'peak_start0','peak_stop_hg38':'peak_end0'})
    require((source.peak_start0 < source.peak_end0).all(), 'Invalid source peak interval.')
    manifest = pd.read_csv(MANIFEST, sep='\t', usecols=[*GEOMETRY_META,'key','chr','pos_hg38','ref','alt','assembly'],
                           dtype=str, keep_default_na=False)
    excluded = pd.read_csv(EXCLUDED, sep='\t', usecols=GEOMETRY_META, dtype=str, keep_default_na=False)
    for table in (manifest, excluded):
        integer(table, ['heldout_fold','length_bp','peak_start0','peak_end0'])
        require(table.heldout_fold.isin([1,2,3,4]).all(), 'Geometry contains forbidden fold.')
        require(table.length_bp.eq(1048576).all(), 'Geometry context changed.')
        require(not table.variant_id.duplicated().any(), 'Repeated geometry variant.')
    integer(manifest, ['pos_hg38'])
    require(len(manifest) == 25477 and len(excluded) == 271, 'Fixed geometry population changed.')
    require(manifest.assembly.eq('GRCh38').all(), 'Manifest assembly changed.')
    require(manifest.coordinate_state.eq('geometry_agrees').all()
            and manifest.target_readout_state.eq('whole_peak_contained').all(), 'Invalid admitted geometry.')
    require(excluded.target_readout_state.ne('whole_peak_contained').all(), 'Excluded whole peak.')
    all_geometry = pd.concat([manifest[GEOMETRY_META], excluded], ignore_index=True)
    require(not all_geometry.variant_id.duplicated().any(), 'Admitted/excluded overlap.')
    require(set(source.lead_variant_id) == set(all_geometry.variant_id), 'Source and entire geometry key sets differ.')
    full = all_geometry.set_index('variant_id').loc[source.lead_variant_id].reset_index(drop=True)
    for column in ('peak_id','heldout_fold','peak_start0','peak_end0'):
        require(np.array_equal(source[column].to_numpy(),full[column].to_numpy()), 'Source/geometry mismatch: '+column)
    require(np.array_equal(source['chr'].str.removeprefix('chr').to_numpy(),
                           full.chromosome.str.removeprefix('chr').to_numpy()), 'Source/geometry chromosome mismatch.')
    wanted = source.set_index('key').loc[manifest.key].reset_index()
    for column in ('key','lead_variant_id','chr','pos_hg38','ref','alt','peak_id','peak_start0','peak_end0','heldout_fold'):
        expected = manifest['variant_id'] if column == 'lead_variant_id' else manifest[column]
        require(np.array_equal(wanted[column].to_numpy(),expected.to_numpy()), 'Manifest/source mismatch: '+column)
    # Recheck immutable files before permitted-effect access.
    for name, path in paths.items():
        require(sha(path) == HASHES[name], 'Input changed during metadata validation: '+name)
    out = Path(os.environ['CODEX_REC_OUTPUT'])/('target_readout_effect_export_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    def save(name, data):
        (out/(name+'.json')).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    save('metadata_admission', {'status':'all_fixed_metadata_checks_passed','input_sha256':HASHES,
                                'source_rows':25748,'manifest_rows':25477,'geometry_exclusions':271,
                                'beta_parsed':False,'p_or_q_values_imported':False})
    # Only now import beta strings, from an already fold-0-free source. Parse
    # numeric beta only for the complete fixed geometry manifest, in its order.
    effects = pd.read_csv(SOURCE, sep='\t', usecols=['lead_variant_id','heldout_fold','beta_alt'],
                          dtype=str, keep_default_na=False)
    require(effects.lead_variant_id.equals(source.lead_variant_id), 'Effect pass source order changed.')
    integer(effects, ['heldout_fold'])
    require(effects.heldout_fold.equals(source.heldout_fold), 'Effect pass fold membership changed.')
    permitted = effects.set_index('lead_variant_id').loc[wanted.lead_variant_id,'beta_alt']
    beta = pd.to_numeric(permitted, errors='raise').to_numpy(float)
    require(np.isfinite(beta).all(), 'Nonfinite admitted source beta.')
    require(sha(SOURCE) == HASHES['source'], 'Source changed during effect read.')
    # Preserve the source decimal text after numeric validation.
    wanted['beta_alt'] = permitted.to_numpy(str)
    output = out/'development_target_effects.tsv.gz'
    wanted[OUTPUT].to_csv(output, sep='\t', index=False, compression={'method':'gzip','mtime':0})
    excluded_keys = set(excluded.variant_id)
    excluded_source = source.loc[source.lead_variant_id.isin(excluded_keys)]
    excluded_source.to_csv(out/'geometry_excluded_metadata.tsv', sep='\t', index=False)
    report = {'status':'development_only_effect_export_complete', 'source':str(SOURCE),
              'manifest':str(MANIFEST),'input_sha256':HASHES,'output':str(output),'output_sha256':sha(output),
              'script_sha256':sha(__file__),'columns':OUTPUT,'source_rows':len(source),
              'covered_rows':len(wanted),'geometry_excluded_rows':len(excluded_source),
              'source_fold_counts':{str(k):int(v) for k,v in source.heldout_fold.value_counts().sort_index().items()},
              'covered_fold_counts':{str(k):int(v) for k,v in wanted.heldout_fold.value_counts().sort_index().items()},
              'excluded_fold_counts':{str(k):int(v) for k,v in excluded_source.heldout_fold.value_counts().sort_index().items()},
              'beta_unit':'unchanged source FastQTL ALT-dosage beta',
              'geometry_unit':'GRCh38 BED0 half-open; variant position one-based',
              'fold0_source_opened':False,'p_or_q_values_imported':False,'metrics_computed':False,
              'fits_performed':False,'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}
    save('summary', report)
    print(json.dumps({'status':report['status'],'covered_rows':len(wanted),
                      'geometry_excluded_rows':len(excluded_source),'metrics_computed':False}),flush=True)


if __name__ == '__main__':
    main()
