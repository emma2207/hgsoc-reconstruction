#!/usr/bin/env nextflow

process DEMULTIPLEXING {
	conda "${params.conda}"
	label "Demultiplexing_${params.pool}"
	publishDir 'demultiplexing/', mode: 'copy'

	input:
	tuple val(pool), path(cellSNP)
	path ref_vcf_file

	output:
	tuple val(pool), path("vireo/pool${pool}")

	script:
	"""
	output_location="vireo/pool${pool}"
	mkdir -p \$output_location

	# Adjust the number of samples for pools 9 & 10
	if ${params.dataset} == "hgsoc-new"
	then
		if [ "$pool" -lt 9 ]; then
			n_samples=4
		elif [ "$pool" -lt 10 ]; then
			n_samples=3
		else
			n_samples=2
		fi
	else
		n_samples=4
	fi

	vireo \
		-c $cellSNP \
		-N \$n_samples \
		-o \$output_location \
		-d $ref_vcf_file \
		--randSeed=12

	"""
}
