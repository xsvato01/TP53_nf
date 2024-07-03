include { TRIMMING_1; TRIMMING_2; ALIGN; SORT_INDEX; PILE_UP;
VARSCAN; VARDICT; NORMALIZE_VARIANTS; MERGE_VARIANTS; NORMALIZE_MERGED_VARIANTS; ANNOTATE;
NORMALIZE_VEP; CREATE_TXT; CREATE_FINAL_TABLE; COVERAGE; COVERAGE_STATS; PLOT_SAMPLE_BARPLOTS;
MERGE_BAMS; BAM_READCOUNT; PLOT_INTERACTIVE_BARPLOTS } from "${params.projectDirectory}/modules"


include { subworkflow_bryja } from "${params.projectDirectory}/subworkflow_bryja"

workflow {
rawfastq = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    . splitCsv( header:true )
    . map { row ->
        def meta = [name:row.name, run:row.run, bryja:row.bryja]
        def baseDir = new File("${params.baseDir}")
		def runDir = baseDir.listFiles(new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.endsWith(meta.run)
			}
		})[0] //get the real folderName that has prepended date
        [meta.name, meta, [
            file("${runDir}/raw_fastq/${meta.name}_R1.fastq.gz", checkIfExists: true),
            file("${runDir}/raw_fastq/${meta.name}_R2.fastq.gz", checkIfExists: true),
        ]]
    }
    // . view()

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
norm_merged	= NORMALIZE_MERGED_VARIANTS(merged)


annotated	= ANNOTATE(norm_merged)
annot_norm	= NORMALIZE_VEP(annotated)
txt		= CREATE_TXT(annot_norm)
runName_samplePath = txt.map({return [it[1].run, it[2]]})
CREATE_FINAL_TABLE(runName_samplePath.groupTuple())

// ///vcf bar plots
norm_merged_whole = norm_merged.map{[it[1].run, it[2]]}
    .groupTuple().map{[it[0], [run:it[0]], it[1]]} //twice run name to be consistent with other samples
PLOT_SAMPLE_BARPLOTS(norm_merged.mix(norm_merged_whole))

//interactive bar plots
groupedBams = sortedbam.map{ [it[1].run, it[2], it[3]]}.groupTuple()
mergedBams = MERGE_BAMS(groupedBams)
readCounts = BAM_READCOUNT(mergedBams)
PLOT_INTERACTIVE_BARPLOTS(readCounts)
///coverage
covered		= COVERAGE(sortedbam)
runName_covPath = covered.map({return [it[1].run, it[2]]})
COVERAGE_STATS(runName_covPath.groupTuple())		

subworkflow_bryja(trimmed2, aligned )	
}