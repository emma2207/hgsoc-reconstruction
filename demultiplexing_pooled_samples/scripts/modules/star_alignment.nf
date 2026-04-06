#!/usr/bin/env nextflow

process STAR_ALIGNMENT {
	conda "${params.conda}"
	tag "pool${pool}"
	publishDir 'demultiplexing/', mode: 'copy'

	input:
	tuple val(pool), path(fastq_r1), path(fastq_r2)

	output:
	tuple val(pool), 
		path("star/${params.dataset}/pool${pool}/Aligned.sortedByCoord.out.bam"), 
		path("star/${params.dataset}/pool${pool}/Aligned.sortedByCoord.out.bam.bai"),
		path("star/${params.dataset}/pool${pool}/Solo.out/Gene/filtered/barcodes.tsv")
	

	script:
	def fastq_r1_arg = fastq_r1.collect { fastq -> fastq.toString() }.join(',')
	def fastq_r2_arg = fastq_r2.collect { fastq -> fastq.toString() }.join(',')
	"""
	# Check if STAR index exists
	if [ -f "${params.refGenome}/star/SAindex" ]
	then
		echo "STAR index exists, skipping."
	else
		STAR \
			--runMode genomeGenerate \
			--runThreadN 6 \
			--genomeDir "${params.refGenome}/star" \
			--genomeFastaFiles "${params.refGenome}/fasta/genome.fa" \
			--sjdbGTFfile "${params.refGenome}/genes/genes.gtf"
	fi

	# Create output directory
	output_location="star/${params.dataset}/pool${pool}"
	mkdir -p \$output_location

	echo "FASTQ R1: ${fastq_r1_arg}"
	echo "FASTQ R2: ${fastq_r2_arg}"
	echo "Reference Genome: ${params.refGenome}"

	# Run STAR alignment
	STAR \
		--outSAMtype BAM SortedByCoordinate \
		--soloType CB_UMI_Simple \
		--soloCBwhitelist None \
		--soloCBstart 1 \
		--soloCBlen 16 \
		--soloUMIstart 17 \
		--soloUMIlen 12 \
		--soloBarcodeReadLength 1 \
		--outSAMattributes NH HI AS nM CB UB \
		--genomeDir "${params.refGenome}/star" \
		--runThreadN 6 \
		--readFilesIn ${fastq_r1_arg} ${fastq_r2_arg} \
		--quantMode GeneCounts \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix "\${output_location}/"

	echo "Start indexing..."

	samtools index \${output_location}/Aligned.sortedByCoord.out.bam

	echo "Finished indexing!"
	"""
}
