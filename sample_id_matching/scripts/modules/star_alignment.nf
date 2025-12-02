#!/usr/bin/env nextflow

process ALIGNMENT_WITH_STAR {

	label 'star'
	publishDir 'results/', mode: 'copy'
	errorStrategy 'ignore'
	cache false

	input:
	tuple val(sample_name), path(fastq_r1), path(fastq_r2)

	output:
	tuple val(sample_name), 
		path("star/${params.dataset}/${params.datatype}/${sample_name}/Aligned.sortedByCoord.out.bam"), 
		path("star/${params.dataset}/${params.datatype}/${sample_name}/Aligned.sortedByCoord.out.bam.bai"), emit: aligned_reads
	path("star/${params.dataset}/${params.datatype}/${sample_name}/ReadsPerGene.out.tab")

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
	output_location="star/${params.dataset}/${params.datatype}/${sample_name}"
	mkdir -p \$output_location

	echo "FASTQ R1: $fastq_r1"
	echo "FASTQ R2: $fastq_r2"
	echo "Reference Genome: ${params.refGenome}"

	# Run STAR alignment
	if [[ ( ${params.datatype} == "single-cell" || ${params.datatype} == "single-nucleus" ) && "${params.read_type}" == "paired" ]]; then
		echo "Running STAR in single-cell mode with paired-end reads"
		STAR \
			--outSAMtype BAM SortedByCoordinate \
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
			--quantMode GeneCounts \
			--readFilesCommand gunzip -c \
			--outFileNamePrefix "\${output_location}/"
	elif [ "${params.read_type}" == "paired" ]
	then
		echo "Running STAR in bulk RNA-seq mode with paired-end reads"
		STAR \
			--outSAMtype BAM SortedByCoordinate \
			--genomeDir ${params.refGenome} \
			--runThreadN 6 \
			--readFilesIn $fastq_r1 $fastq_r2 \
			--readFilesCommand gunzip -c \
			--outFileNamePrefix "\${output_location}/" \
			--quantMode GeneCounts
	elif [ "${params.datatype}" == "single-cell" ] || [ "${params.datatype}" == "single-nucleus" ]
	then
		echo "Running STAR in single-cell mode with single-end reads"
		STAR \
			--outSAMtype BAM SortedByCoordinate \
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
			--readFilesIn $fastq_r1 \
			--quantMode GeneCounts \
			--readFilesCommand gunzip -c \
			--outFileNamePrefix "\${output_location}/"
	else
		echo "Running STAR in bulk RNA-seq mode with single-end reads"
		STAR \
			--outSAMtype BAM SortedByCoordinate \
			--genomeDir ${params.refGenome} \
			--runThreadN 6 \
			--readFilesIn $fastq_r1 \
			--readFilesCommand gunzip -c \
			--outFileNamePrefix "\${output_location}/" \
			--quantMode GeneCounts
	fi
	echo "Start indexing..."

	samtools index \${output_location}/Aligned.sortedByCoord.out.bam

	echo "Finished indexing!"
	"""
}
