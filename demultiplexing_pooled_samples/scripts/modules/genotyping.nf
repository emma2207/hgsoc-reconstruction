#!/usr/bin/env nextflow

process GENOTYPING {
	conda "${params.conda}"
	label "Genotyping_${params.pool}"
	publishDir 'demultiplexing/', mode: 'copy'

	input:
	tuple val(pool), path(bam), path(bai), path(barcodes)

	output:
	tuple val(pool), path("cellSNP/pool${pool}")

	script:
	"""
	output_location="cellSNP/pool${pool}"
	mkdir -p \$output_location

	cellsnp-lite \
		-s $bam \
		-b $barcodes \
		-O \$output_location \
		-p 20 \
		--minMAF=0.1 \
		--minCOUNT=100 \
		--gzip

	"""
}
