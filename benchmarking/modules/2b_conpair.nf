#!/usr/bin/env nextflow

process CONPAIR {
    conda "${params.conda}/conpair"
    publishDir "${params.outdir}/conpair", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(bam)
    
    output:
        path("${params.dataset}/*/*_pileup.txt")
        path("${params.dataset}/*/*_concordance.txt")
    
    script:
        def conpairMarkerBed = params.conpair_marker_bed ?: "${params.CONPAIR}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
        def conpairMarkerTxt = params.conpair_marker_txt ?: "${params.CONPAIR}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
        def conpairMinCov = params.conpair_min_cov ?: 3
        def conpairMinMapQual = params.conpair_min_map_qual ?: 0
        def conpairMinBaseQual = params.conpair_min_base_qual ?: 10

        """
        set -euo pipefail

        export CONPAIR_DIR=/projects/\${USER}/software/conpair
        export GATK_JAR=/projects/\${USER}/software/anaconda/envs/conpair/opt/gatk-3.8/GenomeAnalysisTK.jar
        export PYTHONPATH=/projects/\${USER}/software/conpair/modules

        conpair_marker_bed="${conpairMarkerBed}"
        conpair_marker_txt="${conpairMarkerTxt}"
        conpair_min_cov=${conpairMinCov}
        conpair_min_map_qual=${conpairMinMapQual}
        conpair_min_base_qual=${conpairMinBaseQual}
        conpair_normalize_pileup_py="${params.projectDir}/modules/helpers/conpair_normalize_pileup.py"

        if [ ! -s "\$conpair_marker_bed" ]
        then
            echo "ERROR: Conpair BED marker file not found: \$conpair_marker_bed" >&2
            exit 1
        fi

        if [ ! -s "\$conpair_marker_txt" ]
        then
            echo "ERROR: Conpair TXT marker file not found: \$conpair_marker_txt" >&2
            exit 1
        fi

        if [ ! -s "\$conpair_normalize_pileup_py" ]
        then
            echo "ERROR: Conpair pileup normalizer not found: \$conpair_normalize_pileup_py" >&2
            exit 1
        fi

        fasta_contig="\$(grep -m1 '^>' "${params.refGenome}/fasta/genome.fa" | sed 's/^>//' | awk '{print \$1}')"
        marker_contig="\$(awk 'NF {print \$1; exit}' "\$conpair_marker_txt")"
        remove_chr_opt=""

        if [[ "\$fasta_contig" == chr* && "\$marker_contig" != chr* ]]
        then
            remove_chr_opt="--remove_chr_prefix"
        elif [[ "\$fasta_contig" != chr* && "\$marker_contig" == chr* ]]
        then
            echo "ERROR: Reference FASTA contigs do not use 'chr' prefix, but selected Conpair markers do: \$conpair_marker_txt" >&2
            echo "Set --conpair_marker_bed/--conpair_marker_txt to marker files that match your reference naming." >&2
            exit 1
        fi

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk"
        else
            output_location="${params.dataset}/real_data" 
        fi
        mkdir -p \$output_location

        # Convert bam input to bash array and count files
        bams=($bam)
        count=0
        for file in $bam; do
            count=\$((count + 1))
        done

        # Run conpair over all pairs of bams (without repetition)
        for ((i = 0; i < count; i++))
        do  
            bam_1="\${bams[i]}"
            sample_1=\$(basename "\${bam_1}" .bam)
            pileup_1="\${output_location}/\${sample_1}_pileup.txt"
            pileup_1_raw="\${output_location}/\${sample_1}_pileup.raw.txt"

            python ${params.CONPAIR}/scripts/run_gatk_pileup_for_sample.py \\
                    -B \${bam_1} \\
                    -R ${params.refGenome}/fasta/genome.fa \\
                    -M \$conpair_marker_bed \
                    \$remove_chr_opt \
                    -O \${pileup_1_raw}

                # Normalize raw mpileup to Conpair's expected 5-column format with base/quality alignment.
                python "\$conpair_normalize_pileup_py" "\${pileup_1_raw}" "\${pileup_1}"

            for ((j = i; j < count; j++))
            do
                bam_2="\${bams[j]}"
                sample_2=\$(basename "\${bam_2}" .bam)
                pileup_2="\${output_location}/\${sample_2}_pileup.txt"
                pileup_2_raw="\${output_location}/\${sample_2}_pileup.raw.txt"

                output_file="\${output_location}/\${sample_1}_vs_\${sample_2}_concordance.txt"
                
                # Run CONPAIR
                python ${params.CONPAIR}/scripts/run_gatk_pileup_for_sample.py \\
                    -B \${bam_2} \\
                    -R ${params.refGenome}/fasta/genome.fa \\
                    -M \$conpair_marker_bed \
                    \$remove_chr_opt \
                    -O \${pileup_2_raw}

                python "\$conpair_normalize_pileup_py" "\${pileup_2_raw}" "\${pileup_2}"

                pileup_contig_1="\$(awk 'NF {print \$1; exit}' \${pileup_1})"
                pileup_contig_2="\$(awk 'NF {print \$1; exit}' \${pileup_2})"
                if [[ -z "\$pileup_contig_1" || -z "\$pileup_contig_2" ]]
                then
                    echo "ERROR: Empty pileup detected for one or both samples: \${pileup_1}, \${pileup_2}" >&2
                    exit 1
                fi

                python ${params.CONPAIR}/scripts/verify_concordance.py \\
                    -T \${pileup_1} \
                    -N \${pileup_2} \
                    -M \$conpair_marker_txt \
                    -C \$conpair_min_cov \
                    -Q \$conpair_min_map_qual \
                    -B \$conpair_min_base_qual \
                    -O \${output_file}

                # Conpair exits successfully even when no markers pass filters and may skip writing output.
                # Retry once with permissive thresholds to recover concordance for low-depth datasets.
                if [ ! -s "\${output_file}" ]
                then
                    echo "WARN: No concordance output for \${sample_1} vs \${sample_2} at C=\$conpair_min_cov,Q=\$conpair_min_map_qual,B=\$conpair_min_base_qual; retrying with C=1,Q=0,B=0" >&2
                    python ${params.CONPAIR}/scripts/verify_concordance.py \
                        -T \${pileup_1} \
                        -N \${pileup_2} \
                        -M \$conpair_marker_txt \
                        -C 1 \
                        -Q 0 \
                        -B 0 \
                        -O \${output_file}
                fi

                if [ ! -s "\${output_file}" ]
                then
                    echo "ERROR: Conpair did not produce concordance output for \${sample_1} vs \${sample_2} even after fallback thresholds." >&2
                    echo "Check pileup depth and marker/reference compatibility." >&2
                    exit 1
                fi
            done

        done
        """
}