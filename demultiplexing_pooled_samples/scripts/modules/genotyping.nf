#!/usr/bin/env nextflow

process GENOTYPING {
	conda "${params.conda}"
	tag "pool${pool}"
	publishDir 'demultiplexing/', mode: 'copy'

	input:
	tuple val(pool), path(bam), path(bai), path(barcodes)

	output:
	tuple val(pool), path("cellSNP/pool${pool}")

	script:
	"""
	output_location="cellSNP/pool${pool}"
	mkdir -p \$output_location

	# Remove -1 suffix from barcodes if present (e.g. from 10x v3 chemistry)
	if zcat $barcodes | head -n 1 | grep -qE -- '-1\$'; then
		echo "Removing -1 suffix from barcodes in $barcodes"
		zcat $barcodes | sed 's/-1\$//' | gzip > \$output_location/barcodes.tsv.gz
		stripped_barcodes="\$output_location/barcodes.tsv.gz"
	else
		echo "No -1 suffix found in barcodes, using original file: $barcodes"
		stripped_barcodes=$barcodes
	fi

	cellsnp-lite \
		-s $bam \
		-b \$stripped_barcodes \
		-O \$output_location \
		-p 20 \
		--minMAF=0.1 \
		--minCOUNT=20 \
		--gzip

	"""
}
