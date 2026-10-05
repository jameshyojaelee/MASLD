#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# One explicitly selected library; no array. Original BAM and named counts remain read-only.
set -euo pipefail
: "${SOURCE_BAM:?}" "${SAMPLE_ID:?}" "${LIBRARY_LAYOUT:?}" "${ANNOTATION_ROOT:?}" "${COUNT_OUTPUT_ROOT:?}"
task_root="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas}"
binary_root=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
"$binary_root/python" -B - "$ANNOTATION_ROOT/resolution_summary.json" <<'PY'
import json,sys
s=json.load(open(sys.argv[1]))
assert s.get('annotation_resolution_complete') is True and s.get('n_unresolved_annotations') == 0, 'Full annotation correspondence unresolved'
PY
test -s "$ANNOTATION_ROOT/grouped_counting_annotation.gtf.gz"
test -s "$SOURCE_BAM"
test ! -e "$COUNT_OUTPUT_ROOT"
[[ "$LIBRARY_LAYOUT" == SINGLE || "$LIBRARY_LAYOUT" == PAIRED ]]
scratch_parent="${SLURM_TMPDIR:-${COUNT_SCRATCH_ROOT:?COUNT_SCRATCH_ROOT required without SLURM_TMPDIR}}"
task_scratch=$(mktemp -d "$scratch_parent/resource-count-${SAMPLE_ID}-XXXXXXXX")
flags=(-s 2 -M -R BAM)
if [[ "$LIBRARY_LAYOUT" == PAIRED ]]; then flags+=(-p --countReadPairs -B); fi
assigned="$task_scratch/$(basename "$SOURCE_BAM").featureCounts.bam"
"$binary_root/python" -B - "$binary_root" "$ANNOTATION_ROOT" "$SOURCE_BAM" "$task_scratch" \
  "$SAMPLE_ID" "$LIBRARY_LAYOUT" "${SLURM_CPUS_PER_TASK:-4}" "${flags[@]}" <<'PY'
import hashlib,json,pathlib,platform,subprocess,sys
binary,annotation,source,scratch,sample,layout,cpus,*flags=sys.argv[1:]
source=pathlib.Path(source).resolve(strict=True); scratch=pathlib.Path(scratch)
assert source.name.split('.',1)[0]==sample, 'Source BAM run ID differs from SAMPLE_ID'
summary_path=pathlib.Path(annotation)/'resolution_summary.json'
summary=json.loads(summary_path.read_text())
gtf=pathlib.Path(summary['grouped_counting_annotation']).resolve(strict=True)
correspondence=summary_path.parent/'annotation_correspondence.tsv'
def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1048576),b''): digest.update(block)
    return digest.hexdigest()
assert sha(gtf)==summary['grouped_counting_annotation_sha256']
assert sha(correspondence)==summary['annotation_correspondence_sha256']
def stat(path):
    info=path.stat()
    return {'path':str(path.resolve()),'size_bytes':info.st_size,'mtime_ns':str(info.st_mtime_ns)}
def version(command):
    result=subprocess.run(command,capture_output=True,text=True)
    return (result.stdout+result.stderr).strip()
assigned=scratch/(source.name+'.featureCounts.bam'); sorted_bam=scratch/'name-sorted.bam'
command=[binary+'/featureCounts','-T',cpus,*flags,'-a',str(gtf),'-o',str(scratch/'alignment_counts_not_for_analysis.tsv'),'--Rpath',str(scratch),str(source)]
sort=[binary+'/samtools','sort','-N','-@','2','-m','1G','-T',str(scratch/'name-sort'),'-o',str(sorted_bam),str(assigned)]
evidence={'sample_id':sample,'layout':layout,'source_bam_before':stat(source),
    'featurecounts_argv':command,'samtools_sort_argv':sort,
    'annotation_summary_sha256':sha(summary_path),'grouped_counting_annotation_sha256':sha(gtf),
    'annotation_correspondence_sha256':sha(correspondence),
    'versions':{'python':platform.python_version(),'featureCounts':version([binary+'/featureCounts','-v']),
                'samtools':subprocess.check_output([binary+'/samtools','--version'],text=True).splitlines()[0]}}
completed=subprocess.run(command)
evidence['featurecounts_returncode']=completed.returncode
evidence['source_bam_after']=stat(source)
if completed.returncode==0 and evidence['source_bam_before']==evidence['source_bam_after']:
    evidence['samtools_sort_returncode']=subprocess.run(sort).returncode
    if evidence['samtools_sort_returncode']==0: evidence['assignment_bam']=stat(sorted_bam)
evidence['status']='complete' if evidence.get('samtools_sort_returncode')==0 else 'failed'
with (scratch/'featurecounts_execution.json').open('x') as handle:
    json.dump(evidence,handle,indent=2);handle.write('\n')
if evidence['status']!='complete': raise SystemExit('Assignment execution failed or source BAM changed')
PY
"$binary_root/python" -B "$task_root/scripts/analysis/alt_contig_recount/count_verified_fragments.py" \
  --annotation-summary "$ANNOTATION_ROOT/resolution_summary.json" \
  --assignment-bam "$task_scratch/name-sorted.bam" --sample-id "$SAMPLE_ID" \
  --execution-evidence "$task_scratch/featurecounts_execution.json" \
  --layout "$LIBRARY_LAYOUT" --output "$COUNT_OUTPUT_ROOT"
# Scratch is retained for inspection; scheduler-owned scratch follows its normal lifecycle.
