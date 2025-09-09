#!/usr/bin/env nextflow

process ALIGNMENT_WITH_STAR {

	label 'star'

	input:
	tuple path(fastq_r1), path(fastq_r2)

	output:
	path "star/*"

	stub:
	"""
	echo "Aligning reads for sample ${params.sample} in dataset ${params.dataset}"
	"""

	script:
	"""
	mkdir -p "${params.outputDir}/star/${params.dataset}/${params.sample}"
	output_location="${params.outputDir}/star/${params.dataset}/${params.sample}"

	# Run STAR alignment
	STAR \
		--soloType CB_UMI_Simple \
		--genomeDir ${params.refGenome} \
		--runThreadN 6 \
		--readFilesIn ${fastq_R1} ${fastq_R2} \
		--outFileNamePrefix ${output_location}/ \
		--readFilesCommand gunzip -c \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts

	samtools index ${output_location}/Aligned.sortedByCoord.out.bam
	"""
}
