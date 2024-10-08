include { ALIGN_ANIMALS; FILTER_HUMAN; SORT_INDEX; ZIPFILES; FILESENDER  } from "${params.projectDirectory}/modules"

workflow subworkflow_bryja {
    take:
	trimmed_fastq
	bam

	main:
    animals = Channel.of( 'mouse', 'zebrafish')
    fastqAnimals = trimmed_fastq.filter{it[1].bryja == "TRUE"}.combine(animals).view()
    alignedAnimals = ALIGN_ANIMALS(fastqAnimals)
    patoHuman = bam.filter{it[1].bryja == "TRUE"}
    filteredHuman = FILTER_HUMAN(patoHuman)
    sortedToSendFiles = SORT_INDEX(filteredHuman.mix(alignedAnimals)).map({return [it[2], it[3]]})
    zippedFiles = ZIPFILES(sortedToSendFiles.collect())
    FILESENDER(zippedFiles)
}