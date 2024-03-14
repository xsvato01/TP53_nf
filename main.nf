process COLLECT_BASECALLED {
	tag "COLLECT_BASECALLED on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "small_process"
	label "s_cpu"
	label "s_mem"

	input:
	val(sample)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	"""
	echo COLLECT_BASECALLED $sample.name
	cp  /mnt/shared/MedGen/sequencing_results/primary_data/*${sample.run}/raw_fastq/${sample.name}*R{1,2}* ./
	"""
} 

process TRIMMING_1 {
	tag "trimming 1 on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "s_cpu"
	label "s_mem"
	//publishDir  "${params.outDirectory}/${sample.run}/${sample.name}/trimmed/", mode:'copy'
	
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	"""
	echo TRIMMING_1 $sample.name

	source activate cutadapt
	cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
		-o ${sample.name}.trimmed1.R1.fastq.gz -p ${sample.name}.trimmed1.R2.fastq.gz $reads
	"""
}


process TRIMMING_2 {
	tag "trimming 2 on $sample.name using $task.cpus CPUs and $task.memory memory"
	//publishDir  "${params.outDirectory}/${sample.run}/${sample.name}/trimmed/", mode:'copy'
	label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	""" 
	echo TRIMMING_2 $sample.name
	source activate cutadapt
	cutadapt -a CAAGGGGGACTGTAGATGGG...TAGGATCTGACTGCGGCTCC \
		-A GGAGCCGCAGTCAGATCCTA...CCCATCTACAGTCCCCCTTG \
		-a ACAACGTTCTGGTAAGGACAX -A TGTCCTTACCAGAACGTTGTX --overlap 4 \
		-o ${sample.name}.trimmed2.R1.fastq.gz -p ${sample.name}.trimmed2.R2.fastq.gz $reads
	"""
}

process ALIGN {
	tag "first align on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/${sample.name}/mapped/", mode:'copy'
	label "m_cpu"
	label "l_mem"
	
	input:
	tuple val(sample), path(reads)

	output:
    tuple val(sample), path("${sample.name}.sorted.bam"), path("${sample.name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${sample.name}\\tSM:${sample.name}\\tLB:${sample.name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $sample.name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex} $reads \
	| samtools view -Sb -o - -| samtools sort -o ${sample.name}.sorted.bam
	samtools index ${sample.name}.sorted.bam ${sample.name}.sorted.bai	
	"""
}

process PILE_UP {
	tag "PILE_UP on $sample.name using $task.cpus CPUs and $task.memory memory"
	//publishDir "${params.outDirectory}/${sample.run}/${name}/vcf/", pattern: '*.md.ba*', mode:'copy'
	label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(sample), path(bam), path(bai)

	output:
	tuple val(sample), path("${sample.name}.mpileup")
	
	script:
	"""
	echo PILE_UP $sample.name
	source activate samtools
	samtools mpileup -x -B -Q 25 -d 999999 -L 999999 -F 0.0005 -f ${params.ref}.fa $bam > ${sample.name}.mpileup
	"""
}

process VARSCAN {
	tag "VARSCAN on $sample.name using $task.cpus CPUs and $task.memory memory"
	//publishDir "${params.outDirectory}/${sample.run}/${name}/vcf/", mode:'copy'
    label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(sample), path(mpileup)

	output:
	tuple val(sample), path("${sample.name}.varscan.snv.vcf")
	tuple val(sample), path("${sample.name}.varscan.indel.vcf")
	
	script:
	"""
	echo VARSCAN $sample.name
	source activate varscan
	varscan mpileup2snp $mpileup \
		--strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.000005 \
		--output-vcf | sed 's/Sample1/varscan_SNV/g' > ${sample.name}.varscan.snv.vcf

	varscan mpileup2indel $mpileup \
		--strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.000005 \
		--output-vcf | sed 's/Sample1/varscan_INDEL/g' > ${sample.name}.varscan.indel.vcf
	"""
}

process VARDICT {
	tag "VARDICT on $sample.name using $task.cpus CPUs and $task.memory memory"
	// publishDir "${params.outDirectory}/${sample.run}/${sample.name}/vcf/", mode:'copy'
    label "s_cpu"
	label "m_mem"
	debug true
	
	input:
	tuple val(sample), path(bam), path(bai)

	output:
	tuple val(sample), path ("${sample.name}.vardict.vcf")

	script:
	"""
	echo VARDICT $sample.name
	source activate vardict
	export PATH="/opt/conda/envs/samtools/bin/samtools:$PATH"
	
	sleep infinity
	vardict -G ${params.ref}.fa -f 0.0005 -b ${bam} \
		-c 1 -S 2 -E 3 -r 8 -Q 1 -q 25 -P 2 -m 8 \
		${params.varbed} | Rscript --vanilla ${params.teststrandbias} | perl ${params.var2vcf_valid} -f 0.000005 -d 50 -c 5 -p 2 -q 25 -Q 1 -v 8 -m 8 -N vardict - > ${sample.name}.vardict.vcf
	"""
}



process NORMALIZE_VARIANTS {
	tag "Normalizing variants on $sample.name using $task.cpus CPUs and $task.memory memory"
	// publishDir "${params.outDirectory}/${sample.run}/${sample.name}/vcf/", mode:'copy'
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path (varscan_snv)
	tuple val(sample), path (varscan_indel)
	tuple val(sample), path (vardict)

	output:
	tuple val(sample), path ("${sample.name}.varscan.snv.norm.vcf")
	tuple val(sample), path ("${sample.name}.varscan.indel.norm.vcf")
	tuple val(sample), path ("${sample.name}.vardict.norm.vcf")
	
	script:
	"""
	echo NORMALIZE_VARIANTS $sample.name
	source activate bcftools
	bcftools norm -f ${params.ref}.fa $varscan_snv -o ${sample.name}.varscan.snv.norm.vcf
	bcftools norm -f ${params.ref}.fa $varscan_indel -o ${sample.name}.varscan.indel.norm.vcf
	bcftools norm -f ${params.ref}.fa $vardict -o ${sample.name}.vardict.norm.vcf	
	"""
}


process MERGE_VARIANTS {
	tag "Merging variants on $sample.name using $task.cpus CPUs and $task.memory memory"
	//publishDir "${params.outDirectory}/${sample.run}/${name}/vcf/", mode:'copy'
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path(varscan_snv_norm)
	tuple val(sample), path(varscan_indel_norm)
	tuple val(sample), path(vardict_norm)

	output:
	tuple val(sample), path ("${name}.allcallers.merged.vcf")
	
	script:
	"""
	echo MERGE_VARIANTS $sample.name
	source activate java
	java -jar ${params.gatk36} -T CombineVariants --variant:vardict ${vardict_norm} \
		--variant:varscan_SNV ${varscan_snv_norm} \
		--variant:varscan_INDEL ${varscan_indel_norm} \
		-R ${params.ref}.fa -genotypeMergeOptions UNSORTED \
		--disable_auto_index_creation_and_locking_when_reading_rods -o ${sample.name}.allcallers.merged.vcf
	"""
}



process NORMALIZE_MERGED_VARIANTS {
	tag "Normalizing merged variants on $sample.name using $task.cpus CPUs and $task.memory memory"
	//publishDir "${params.outDirectory}/${sample.run}/${name}/vcf/", mode:'copy'
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path (merged_vcf)

	output:
	tuple val(sample), path ("${sample.name}.allmerged.norm.vcf")
	
	script:
	"""
	echo NORMALIZE_MERGED_VARIANTS $sample.name
	source activate bcftools
	bcftools norm -f ${params.ref}.fa $merged_vcf -o ${sample.name}.allmerged.norm.vcf
	"""
}


process ANNOTATE {
	tag "Annotating variants on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/${sample.name}/annotate/", mode:'copy'
    label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(sample), path (merged_normed_vcf)

	output:
	tuple val(sample), path ("*allmerged.norm.annot.vcf")
	
	script:
	"""
	echo ANNOTATE $sample.name
	source activate vep
	vep -i $merged_normed_vcf --cache --cache_version 87 --dir_cache $params.vep \
	--fasta ${params.ref}.fa --merged --offline --vcf --hgvs --sift b --polyphen b -o ${sample.name}.allmerged.norm.annot.vcf
	"""
}



process NORMALIZE_VEP {
	tag "NORMALIZE_VEP on $sample.name using $task.cpus CPUs and $task.memory memory"
	//publishDir "${params.outDirectory}/${sample.run}/${name}/annotate/", mode:'copy'
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path(annotated)

	output:
	tuple val(sample), path("${sample.name}.allmerged.norm.annot.norm.vcf")
	
	script:
	"""
	echo NORMALIZE_VEP $sample.name
	source activate biopet
	biopet tool VepNormalizer -I $annotated -O ${sample.name}.allmerged.norm.annot.norm.vcf -m explode  
	"""
}


process CREATE_TXT {
	tag "CREATE_TXT on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/annotate/", mode:'copy'
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path(annotated_normed)

	output:
	tuple val(sample), path("${sample.name}.norm.merged.annot.normVEP.txt")
	
	script:
	"""
	echo CREATE_TXT $sample.name
	source activate py36
	python ${params.vcf_simplify} SimplifyVCF -toType table -inVCF $annotated_normed -out ${sample.name}.norm.merged.annot.normVEP.txt  
	"""
}


process CREATE_FINAL_TABLE {
	tag "CREATE_FINAL_TABLE on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/annotate/", mode:'copy'
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path(annotated_normed)
	
	script:
	"""
	echo CREATE_FINAL_TABLE $sample.name
	source activate erko
	Rscript --vanilla ${params.create_table} ${launchDir}/annotate $sample.run  
	"""
}

process CREATE_MERGED_TABLE {
	tag "CREATE_MERGED_TABLE using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/annotate/", mode:'copy', pattern: "*sample.merged.anot.txt"
    publishDir "${params.outDirectory}/${sample.run}/annotate/", mode:'copy', pattern: "*allsamples.merged.anot.txt"
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(sample), path(all_annotated_normed)

    output:
    path "*sample.merged.anot.txt"
    path "*allsamples.merged.anot.txt"

	
	script:
	"""
	echo CREATE_MERGED_TABLE $sample.name
	source activate erko
	Rscript --vanilla ${params.create_table_merged} . $run  
	"""
}


process COVERAGE {
	tag "Creating coverage on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'
    label "s_cpu"
	label "l_mem"
	
	input:
	tuple val(sample), path(bam)

	output:
	tuple val(sample), path("*.txt")
	
	script:
	"""
	echo COVERAGE $sample.name
	source activate bedtools
	bedtools coverage -abam ${params.covbed} -b $bam -d > ${sample.name}.PBcoverage.txt   
	"""
}

process COVERAGE_STATS {
	tag "Creating coverage stats on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'
    label "s_cpu"
	label "l_mem"
	
	input:
	tuple val(sample), path(whatever_navaznost)
	
	script:
	"""
	echo COVERAGE_STATS $sample.name
	source activate erko
	Rscript --vanilla ${params.coverstat} ${launchDir}/coverage $sample.run  
	"""
}


process MULTIQC {
	tag "first QC on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/${sample.name}/coverage/", mode:'copy'
    label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(sample), path(bam)

	output:
	path "report.html"

	script:
	"""
	samtools flagstat $bam > ${sample.name}.flagstat
	samtools stats $bam > ${sample.name}.samstats
    picard BedToIntervalList -I ${params.covbedpicard} -O ${sample.name}.interval_list -SD ${params.ref}.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${sample.name}.interval_list -TARGET_INTERVALS ${sample.name}.interval_list -R ${params.ref}.fa -O ${sample.name}.aln_metrics
	multiqc . -n report.html
	"""

}


workflow {
runlist = channel.fromList(params.samples)
rawfastq = COLLECT_BASECALLED(runlist)

trimmed1	= TRIMMING_1(rawfastq)
trimmed2	= TRIMMING_2(trimmed1)	
sortedbam	= ALIGN(trimmed2)
pileup		= PILE_UP(sortedbam)
varscanned	= VARSCAN(pileup)
vardicted	= VARDICT(sortedbam)
normalized	= NORMALIZE_VARIANTS(varscanned,vardicted)
merged		= MERGE_VARIANTS(normalized)
norm_merged	= NORMALIZE_MERGED_VARIANTS(merged)
// annotated	= ANNOTATE(norm_merged)
// annot_norm	= NORMALIZE_VEP(annotated)
// txt		= CREATE_TXT(annot_norm)
// final_table	= CREATE_FINAL_TABLE(txt)

// covered		= COVERAGE(sortedbam)
// COVERAGE_STATS(covered)					
// MULTIQC(sortedbam)	
}