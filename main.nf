include { REFORMAT_SAMPLE; COLLECT_BASECALLED; TRIMMING_1; TRIMMING_2; ALIGN; SORT_INDEX; PILE_UP;
VARSCAN; VARDICT; NORMALIZE_VARIANTS; MERGE_VARIANTS; NORMALIZE_MERGED_VARIANTS; ANNOTATE;
NORMALIZE_VEP; CREATE_TXT; CREATE_FINAL_TABLE; COVERAGE; COVERAGE_STATS } from "${params.projectDirectory}/modules"


include { subworkflow_bryja } from "${params.projectDirectory}/subworkflow_bryja"

workflow {
runlist = channel.fromList(params.samples)
reformatedInput	= REFORMAT_SAMPLE(runlist)
rawfastq = COLLECT_BASECALLED(reformatedInput)
trimmed1	= TRIMMING_1(rawfastq)
trimmed2	= TRIMMING_2(trimmed1)	
aligned	= ALIGN(trimmed2)
sortedbam	= SORT_INDEX(aligned)
(pileup, _)		= PILE_UP(sortedbam)
varscanned	= VARSCAN(pileup)
vardicted	= VARDICT(sortedbam)
normalized	= NORMALIZE_VARIANTS(varscanned,vardicted)
merged		= MERGE_VARIANTS(normalized)
norm_merged	= NORMALIZE_MERGED_VARIANTS(merged)
annotated	= ANNOTATE(norm_merged)
annot_norm	= NORMALIZE_VEP(annotated)
txt		= CREATE_TXT(annot_norm)
runName_samplePath = txt.map({return [it[1].run, it[2]]})
CREATE_FINAL_TABLE(runName_samplePath.groupTuple())

covered		= COVERAGE(sortedbam)
runName_covPath = covered.map({return [it[1].run, it[2]]})

runName_covPath.groupTuple().view()
COVERAGE_STATS(runName_covPath.groupTuple())		

subworkflow_bryja(trimmed2, aligned )	
}