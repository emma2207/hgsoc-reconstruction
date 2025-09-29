#!/usr/bin/env nextflow

process ALIGNMENT_WITH_STAR {

	label 'star'

	input:
	tuple path(fastq_r1), path(fastq_r2)

	output:
	tuple path("Aligned.sortedByCoord.out.bam"), path("Aligned.sortedByCoord.out.bam.bai")

	script:
	"""

	# Check if STAR index exists
	if [ -f "${params.refGenome}/SAindex" ]
	then
		echo "STAR index exists, skipping."
	else
		STAR \
			--runMode genomeGenerate \
			--runThreadN 20 \
			--genomeDir ${params.refGenome} \
			--genomeFastaFiles "${params.refGenome}/../fasta/genome.fa" \
			--sjdbGTFfile "${params.refGenome}/../genes/genes.gtf"
	fi

	# Create output directory
	mkdir -p "${params.outputDir}/star/${params.dataset}/${params.sample}"
	output_location="${params.outputDir}/star/${params.dataset}/${params.sample}"

	echo "FASTQ R1: $fastq_r1"
	echo "FASTQ R2: $fastq_r2"
	echo "Reference Genome: ${params.refGenome}"

	# Run STAR alignment
	STAR \
		--soloType CB_UMI_Simple \
		--soloCBwhitelist None \
		--soloCBstart 1 \
		--soloCBlen 16 \
        --soloUMIstart 17 \
		--soloUMIlen 12 \
        --soloBarcodeReadLength 0 \
		--outSAMattributes NH HI AS nM CB UB \
		--genomeDir ${params.refGenome} \
		--runThreadN 6 \
		--readFilesIn $fastq_r1 $fastq_r2 \
		--readFilesCommand gunzip -c \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts

	echo "Start indexing..."

	samtools index Aligned.sortedByCoord.out.bam

	cp -r ./* \${output_location}/

	echo "Finished indexing!"
	"""

	stub:
	"""
	echo "Aligning reads for sample ${params.sample} in dataset ${params.dataset}"
	"""
}
