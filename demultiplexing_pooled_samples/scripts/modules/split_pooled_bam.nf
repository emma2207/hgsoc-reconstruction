#!/usr/bin/env nextflow

process SPLIT_POOLED_BAM {
	conda "${params.conda}"
	tag "pool${pool}"
	publishDir 'demultiplexing/', mode: 'copy'

	input:
	tuple val(pool), path(bam), path(bai), path(vireo_pool)

	output:
	tuple val(pool), path("split_bam/pool${pool}")

	script:
	"""
	output_location="split_bam/pool${pool}"
	mkdir -p \$output_location

	donor_file="${vireo_pool}/donor_ids.tsv"
	if [ ! -f "\$donor_file" ]; then
		echo "ERROR: Could not find donor assignments at \$donor_file"
		exit 1
	fi

	awk -F '\t' 'NR>1 && \$2 ~ /^donor/ {print \$2}' "\$donor_file" | sort -u > "\$output_location/donors.txt"

	while IFS= read -r donor
	do
		[ -z "\$donor" ] && continue

		barcode_file="\$output_location/\${donor}_barcodes.txt"
		awk -F '\t' -v donor="\$donor" 'NR>1 && \$2==donor {print \$1}' "\$donor_file" > "\$barcode_file"

		if [ ! -s "\$barcode_file" ]; then
			echo "Skipping \$donor because no barcodes were found"
			continue
		fi

		sample_bam="\$output_location/\${donor}.bam"
		samtools view -D CB,"\$barcode_file" -b "$bam" -o "\$sample_bam"
		samtools index "\$sample_bam"
	done < "\$output_location/donors.txt"
	"""
}
