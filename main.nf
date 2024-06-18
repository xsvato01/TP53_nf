include { REFORMAT_SAMPLE; COLLECT_BASECALLED; TRIMMING_1; TRIMMING_2; ALIGN; SORT_INDEX; PILE_UP;
VARSCAN; VARDICT; NORMALIZE_VARIANTS; MERGE_VARIANTS; NORMALIZE_MERGED_VARIANTS; ANNOTATE;
NORMALIZE_VEP; CREATE_TXT; CREATE_FINAL_TABLE; COVERAGE; COVERAGE_STATS; PLOT_SAMPLE_BARPLOTS } from "${params.projectDirectory}/modules"


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
(varscannedSNV, varscannedINDEL)	= VARSCAN(pileup)
vardicted	= VARDICT(sortedbam)

mixVCFs = vardicted.join(varscannedSNV, by: [0,1]).join(varscannedINDEL, by: [0,1])
normalized	= NORMALIZE_VARIANTS(mixVCFs)
merged		= MERGE_VARIANTS(normalized)
norm_merged	= NORMALIZE_MERGED_VARIANTS(merged)//.view()
norm_merged_whole = norm_merged.map{[it[1].run, it[2]]}
    .groupTuple().map{[it[0], [run:it[0]], it[1]]} //twice run name to be consistent with other samples
norm_merged.mix(norm_merged_whole)
PLOT_SAMPLE_BARPLOTS(norm_merged.mix(norm_merged_whole))
// annotated	= ANNOTATE(norm_merged)
// annot_norm	= NORMALIZE_VEP(annotated)
// txt		= CREATE_TXT(annot_norm)
// runName_samplePath = txt.map({return [it[1].run, it[2]]})

// CREATE_FINAL_TABLE(runName_samplePath.groupTuple())

// covered		= COVERAGE(sortedbam)
// runName_covPath = covered.map({return [it[1].run, it[2]]})
// COVERAGE_STATS(runName_covPath.groupTuple())		

// subworkflow_bryja(trimmed2, aligned )	
}