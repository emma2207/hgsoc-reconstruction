#!/usr/bin/env nextflow

process STAR_ALIGNMENT {
	conda "${params.conda}"
	label "STAR_${params.pool}"
	publishDir 'demultiplexing/', mode: 'copy'

	input:
	tuple val(pool), path(fastq_r1), path(fastq_r2)

	output:
	tuple val(pool), 
		path("star/${params.dataset}/pool${pool}/Aligned.sortedByCoord.out.bam"), 
		path("star/${params.dataset}/pool${pool}/Aligned.sortedByCoord.out.bam.bai")

	script:
	"""
	# Check if STAR index exists
	if [ -f "${params.refGenome}/SAindex" ]
	then
		echo "STAR index exists, skipping."
	else
		STAR \
			--runMode genomeGenerate \
			--runThreadN 6 \
			--genomeDir ${params.refGenome} \
			--genomeFastaFiles "${params.refGenome}/../fasta/genome.fa" \
			--sjdbGTFfile "${params.refGenome}/../genes/genes.gtf"
	fi

	# Create output directory
	output_location="star/${params.dataset}/pool${pool}"
	mkdir -p \$output_location

	echo "FASTQ R1: $fastq_r1"
	echo "FASTQ R2: $fastq_r2"
	echo "Reference Genome: ${params.refGenome}"

	whitelist="${params.whitelist}"
	if [ -f "\$whitelist" ]; then
		echo "Using barcode whitelist: \$whitelist"
	elif [ -f "${params.refGenome}/../\$whitelist" ]; then
		whitelist="${params.refGenome}/../\$whitelist"
		echo "Using barcode whitelist: \$whitelist"
	else
		echo "WARNING: barcode whitelist not found ('\$whitelist'); falling back to '--soloCBwhitelist None'"
		whitelist="None"
	fi

	# Run STAR alignment
	STAR \
		--outSAMtype BAM SortedByCoordinate \
		--soloType CB_UMI_Simple \
		--soloCBwhitelist \$whitelist \
		--soloCBstart 1 \
		--soloCBlen 16 \
		--soloUMIstart 17 \
		--soloUMIlen 12 \
		--soloBarcodeReadLength 0 \
		--outSAMattributes NH HI AS nM CB UB \
		--genomeDir ${params.refGenome} \
		--runThreadN 6 \
		--readFilesIn $fastq_r1 $fastq_r2 \
		--quantMode GeneCounts \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix "\${output_location}/"

	echo "Start indexing..."

	samtools index \${output_location}/Aligned.sortedByCoord.out.bam

	echo "Finished indexing!"
	"""
}
