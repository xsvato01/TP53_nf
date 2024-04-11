process REFORMAT_SAMPLE {
	tag "Reformating $sample.name using $task.cpus CPUs $task.memory"
	label "s_cpu"
	label "xxs_mem"

	input:
	val sample

	output:
	tuple val(sample.name), val(sample)

	""" """ //this is not an error
}

process COLLECT_BASECALLED {
	tag "COLLECT_BASECALLED on $name using $task.cpus CPUs and $task.memory memory"
	label "s_cpu"
	label "xxs_mem"

	input:
	tuple val(name), val(sample)

	output:
	tuple val(name), val(sample), path("*.fastq.gz")

	script:
	"""
	echo COLLECT_BASECALLED $name
	cp  /mnt/share/710000-CEITEC/713000-cmm/713003-pospisilova/base/sequencing_results/primary_data/*${sample.run}/raw_fastq/${name}*R{1,2}* ./
	"""
} 

process TRIMMING_1 {
	tag "trimming 1 on $name using $task.cpus CPUs and $task.memory memory"
	label "s_cpu"
	label "xxs_mem"
	//publishDir  "${params.outDirectory}/${sample.run}/${name}/trimmed/", mode:'copy'
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
	tuple val(name), val(sample), path("*.fastq.gz")

	script:
	"""
	echo TRIMMING_1 $name

	source activate cutadapt
	cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
		-o ${name}.trimmed1.R1.fastq.gz -p ${name}.trimmed1.R2.fastq.gz $reads
	"""
}


process TRIMMING_2 {
	tag "trimming 2 on $name using $task.cpus CPUs and $task.memory memory"
	//publishDir  "${params.outDirectory}/${sample.run}/${name}/trimmed/", mode:'copy'
	label "s_cpu"
	label "xxs_mem"
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
	tuple val(name), val(sample), path("*.fastq.gz")

	script:
	""" 
	echo TRIMMING_2 $name
	source activate cutadapt
	cutadapt -a CAAGGGGGACTGTAGATGGG...TAGGATCTGACTGCGGCTCC \
		-A GGAGCCGCAGTCAGATCCTA...CCCATCTACAGTCCCCCTTG \
		-a ACAACGTTCTGGTAAGGACAX -A TGTCCTTACCAGAACGTTGTX --overlap 4 \
		-o ${name}.trimmed2.R1.fastq.gz -p ${name}.trimmed2.R2.fastq.gz $reads
	"""
}

process ALIGN {
	tag "ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	// publishDir "${params.outDirectory}/${sample.run}/${name}/mapped/", mode:'copy'
	label "m_cpu"
	label "l_mem"
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
    tuple val(name), val(sample), path("${name}.bam")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN $name
	source activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex} $reads > ${name}.sam
	samtools view -Sb ${name}.sam -o ${name}.bam
	"""
}

process FILTER_HUMAN {
	tag "FILTER_HUMAN on $name using $task.cpus CPUs and $task.memory memory"
	// publishDir "${params.outDirectory}/${sample.run}/${name}/mapped/", mode:'copy'
	label "xxs_mem"
	label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam)

	output:
	tuple val("${name}.human"), val(sample), path("${name}.human.bam")


	script:
	"""
	echo FILTER_HUMAN $name
	source activate bedtools
	bedtools intersect -v -abam ${bam} -b $params.varbed > ${name}.human.bam
	"""
}

process ALIGN_ANIMALS {
	tag "ALIGN_ANIMALS on $name to $reference using $task.cpus CPUs and $task.memory memory"
	// publishDir "${params.outDirectory}/${sample.run}/${name}/mapped/", mode:'copy'
	label "m_cpu"
	memory "${reference == 'zebrafish' ? '96 GB' : '48 GB'}"
	
	input:
	tuple val(name), val(sample), path(reads), val(reference)

	output:
    tuple val("${name}.${reference}"), val(sample), path("${name}.${reference}.bam")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo ALIGN_ANIMALS $name to $reference
	source activate bwa
	if [ ${reference} = "mouse" ]; then
		bwa mem -R ${rg} -t $task.cpus ${params.bwa_index_mouse} $reads | samtools view -b -F 4 -o ${name}.${reference}.bam
	elif [ ${reference} = "zebrafish" ]; then
		bwa mem -R ${rg} -t $task.cpus ${params.bwa_index_zebrafish} $reads | samtools view -b -F 4 -o ${name}.${reference}.bam
	fi
	"""
}

process SORT_INDEX {
	tag "Sort index on $name using $task.cpus CPUs and $task.memory memory"
	label "m_mem"
	label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam)

	output:
	tuple val(name), val(sample), path("${name}.sorted.bam"), path("${name}.sorted.bai")


	script:
	"""
	echo SORT_INDEX $name
	source activate samtools
	samtools sort $bam -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai
	"""
}

process PILE_UP {
	tag "PILE_UP on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy', pattern: "*.ba*"
	label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.mpileup")
	tuple path("$bam"), path("$bai")
	
	script:
	"""
	echo PILE_UP $name
	source activate samtools
	samtools mpileup -x -B -Q 25 -d 999999 -L 999999 -F 0.0005 -f ${params.ref}.fa $bam > ${name}.mpileup
	"""
}

process VARSCAN {
	tag "VARSCAN on $name using $task.cpus CPUs and $task.memory memory"
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(name), val(sample), path(mpileup)

	output:
	tuple val(name), val(sample), path("${name}.varscan.snv.vcf")
	tuple val(name), val(sample), path("${name}.varscan.indel.vcf")
	
	script:
	"""
	echo VARSCAN $name
	source activate varscan
	varscan mpileup2snp $mpileup \
		--strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.000005 \
		--output-vcf | sed 's/Sample1/varscan_SNV/g' > ${name}.varscan.snv.vcf

	varscan mpileup2indel $mpileup \
		--strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.000005 \
		--output-vcf | sed 's/Sample1/varscan_INDEL/g' > ${name}.varscan.indel.vcf
	"""
}

process VARDICT {
	tag "VARDICT on $name using $task.cpus CPUs and $task.memory memory"
    label "s_cpu"
	label "xs_mem"
	
	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.vardict.vcf")

	script:
	"""
	echo VARDICT $name
	source activate vardict
	vardict -G ${params.ref}.fa -f 0.0005 -b ${bam} \
		-c 1 -S 2 -E 3 -r 8 -Q 1 -q 25 -P 2 -m 8 \
		${params.varbed} | Rscript --vanilla ${params.teststrandbias} | perl ${params.var2vcf_valid} -f 0.000005 -d 50 -c 5 -p 2 -q 25 -Q 1 -v 8 -m 8 -N vardict - > ${name}.vardict.vcf
	"""
}



process NORMALIZE_VARIANTS {
	tag "Normalizing variants on $name using $task.cpus CPUs and $task.memory memory"
	container "staphb/bcftools:1.10.2"
    label "s_cpu"
	label "xxs_mem"
	
	input:
	tuple val(name), val(sample), path(vardict), path(varscan_snv), path(varscan_indel)

	output:
	tuple val(name), val(sample), path("${name}.varscan.snv.norm.vcf"), path("${name}.varscan.indel.norm.vcf"), path("${name}.vardict.norm.vcf")
	
	script:
	"""
	echo NORMALIZE_VARIANTS $name
	bcftools norm -f ${params.ref}.fa $varscan_snv -o ${name}.varscan.snv.norm.vcf
	bcftools norm -f ${params.ref}.fa $varscan_indel -o ${name}.varscan.indel.norm.vcf
	bcftools norm -f ${params.ref}.fa $vardict -o ${name}.vardict.norm.vcf	
	"""
}


process MERGE_VARIANTS {
	tag "Merging variants on $name using $task.cpus CPUs and $task.memory memory"
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(name), val(sample), path(varscan_snv_norm), path(varscan_indel_norm), path(vardict_norm)

	output:
	tuple val(name), val(sample), path("${name}.allcallers.merged.vcf")
	
	script:
	"""
	echo MERGE_VARIANTS $name
	source activate java
	java -jar ${params.gatk36} -T CombineVariants --variant:vardict ${vardict_norm} \
		--variant:varscan_SNV ${varscan_snv_norm} \
		--variant:varscan_INDEL ${varscan_indel_norm} \
		-R ${params.ref}.fa -genotypeMergeOptions UNSORTED \
		--disable_auto_index_creation_and_locking_when_reading_rods -o ${name}.allcallers.merged.vcf
	"""
}



process NORMALIZE_MERGED_VARIANTS {
	tag "Normalizing merged variants on $name using $task.cpus CPUs and $task.memory memory"
	container "staphb/bcftools:1.10.2"
    label "s_cpu"
	label "xxs_mem"
	
	input:
	tuple val(name), val(sample), path(merged_vcf)

	output:
	tuple val(name), val(sample), path("${name}.allmerged.norm.vcf")
	
	script:
	"""
	echo NORMALIZE_MERGED_VARIANTS $name
	bcftools norm -f ${params.ref}.fa $merged_vcf -o ${name}.allmerged.norm.vcf
	"""
}


process ANNOTATE {
	tag "Annotating variants on $name using $task.cpus CPUs and $task.memory memory"
    label "s_cpu"
	label "s_mem"
	
	input:
	tuple val(name), val(sample), path(merged_normed_vcf)

	output:
	tuple val(name), val(sample), path("${name}.allmerged.norm.annot.vcf")
	
	script:
	"""
	echo ANNOTATE $name
	source activate vep
	vep -i $merged_normed_vcf --cache --cache_version 87 --dir_cache $params.vep \
	--fasta ${params.ref}.fa --merged --offline --vcf --hgvs --sift b --polyphen b -o ${name}.allmerged.norm.annot.vcf
	"""
}



process NORMALIZE_VEP {
	tag "NORMALIZE_VEP on $name using $task.cpus CPUs and $task.memory memory"
    label "s_cpu"
	label "xs_mem"
	
	input:
	tuple val(name), val(sample), path(annotated)

	output:
	tuple val(name), val(sample), path("${name}.allmerged.norm.annot.norm.vcf")
	
	script:
	"""
	echo NORMALIZE_VEP $name
	source activate biopet
	biopet tool VepNormalizer -I $annotated -O ${name}.allmerged.norm.annot.norm.vcf -m explode  
	"""
}


process CREATE_TXT {
	tag "CREATE_TXT on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/tmp/", mode:'copy'

    label "s_cpu"
	label "xxs_mem"
	
	input:
	tuple val(name), val(sample), path(annotated_normed)

	output:
	tuple val(name), val(sample), path("${name}.norm.merged.annot.normVEP.txt")
	
	script:
	"""
	echo CREATE_TXT $name
	source activate py36
	python ${params.vcf_simplify} SimplifyVCF -toType table -inVCF $annotated_normed -out ${name}.norm.merged.annot.normVEP.txt  
	"""
}


process CREATE_FINAL_TABLE {
	tag "CREATE_FINAL_TABLE on $run using $task.cpus CPUs and $task.memory memory"
    publishDir "${params.outDirectory}/${run}/annotate/", mode:'copy'
	label "m_mem"
    label "s_cpu"
	
	input:
	tuple val(run), path(all_annotated_normed)

    output:
    path "*sample.merged.anot.txt"
    path "*allsamples.merged.anot.txt"
	
	script:
	"""
	echo CREATE_FINAL_TABLE $run
	source activate erko
	Rscript --vanilla ${params.create_table} . $run
	"""
}

process COVERAGE {
	tag "COVERAGE on $name using $task.cpus CPUs and $task.memory memory"
	// publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'
    label "s_cpu"
	label "l_mem"
	
	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.PBcoverage.txt")
	
	script:
	"""
	echo COVERAGE $name
	source activate bedtools
	bedtools coverage -abam ${params.covbed} -b $bam -d > ${name}.PBcoverage.txt   
	"""
}

process COVERAGE_STATS {
	tag "COVERAGE_STATS on $run using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${run}/coverage/", mode:'copy'
    label "s_cpu"
	label "xs_mem"
	
	input:
	tuple val(run), path(all_annotated_normed)

	output:
    path "*.perexon_stat.txt"
    path "*_tp53_gene_coverage.txt"
	path "*.allsamples.merged.coverage.txt"
	
	script:
	"""
	echo COVERAGE_STATS $run
	source activate erko
	Rscript --vanilla ${params.coverstat} . $run  
	"""
}


process MULTIQC {
	tag "MULTIQC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/${name}/coverage/", mode:'copy'
    label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(name), val(sample), path(bam)

	output:
	path"report.html"

	script:
	"""
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
    picard BedToIntervalList -I ${params.covbedpicard} -O ${name}.interval_list -SD ${params.ref}.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${name}.interval_list -TARGET_INTERVALS ${name}.interval_list -R ${params.ref}.fa -O ${name}.aln_metrics
	multiqc . -n report.html
	"""
}


process ZIPFILES {
	tag "ZIPFILES $task.cpus CPUs and $task.memory memory"
	// publishDir "${params.outDirectory}/zip/", mode:'copy'
    label "s_cpu"
	label "m_mem"
	
	input:
	pathfiles

	output:
	path"NGS.tar.gz"

	script:
	"""
	echo ZIPFILES
	tar -chf NGS.tar.gz --exclude='.command*' .
	"""
}

process FILESENDER {
	tag "FILESENDER $task.cpus CPUs and $task.memory memory"
    label "s_cpu"
	label "m_mem"
	
	input:
	pathfileToSend

	script:
	"""
	python3 $params.filesender -a 53ee39671b0b915c2f393a75966f0b74683eb0b3b02b6385da604a946689a86d \
	-u 8a7faace1e24c4189f4e3f51d7a8713555a18540@einfra.cesnet.cz -r hana.plesingerova@sci.muni.cz $fileToSend
	"""
}