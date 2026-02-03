#!/usr/bin/env nextflow

process BAMIXCHECKER {
    conda "${params.conda}/bamixchecker"
    publishDir "${params.outdir}/2b_bamixchecker", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(bam)
    
    output:
        path("${params.dataset}/*/BAMixChecker/BAMixChecker_Report.html")
        path("${params.dataset}/*/BAMixChecker/BAMixChecker_Heatmap.pdf")
        path("${params.dataset}/*/BAMixChecker/Total_result.txt")
        path("${params.dataset}/*/BAMixChecker/Matched_samples.txt"), optional: true
        path("${params.dataset}/*/BAMixChecker/Mismatched_samples.txt"), optional: true
    script:
        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}"
        else
            output_location="${params.dataset}/real_data" 
        fi

        # Create config file for BAMixChecker
        echo "GATK=${params.conda}/bamixchecker/bin/gatk" > BAMixChecker.config
        echo "BEDTOOLS=${params.conda}/bamixchecker/bin/bedtools" >> BAMixChecker.config

        # Fix the read groups
        for file in ${bam}
        do  
            base_name=\$(basename "\$file")
            base_name=\${base_name%.bam}

            # Sort reads numerically instead of lexicographically
            # extract header and isolate ##contig lines
            # samtools view -H "\${file}" > "header.txt"
            #grep '^@HD' "header.txt" > "start.txt"
            #grep '^@SQ' "header.txt" > "contigs.txt"
            #grep -e '^@PG' -e '^@CO' "header.txt" > "end.txt"
            
            # produce sortable keys for contigs: numeric chrN -> N, chrX->50000, chrY->50001, chrM/MT->50002, others->100000+name
            #awk '
            #{
            #    if (match(\$0, /SN:[^[:space:]]+/)) {
            #    sn = substr(\$0, RSTART + 3, RLENGTH - 3)
            #    } else {
            #    sn = ""
            #    }
            #    if (sn ~ /^chr[0-9]+\$/) {
            #    n = substr(sn,4) + 0
            #    printf("%05d\\t\\t%s\\n", n, \$0)
            #    } else if (sn == "chrX") {
            #    printf("%05d\\t\\t%s\\n", 50000, \$0)
            #    } else if (sn == "chrY") {
            #    printf("%05d\\t\\t%s\\n", 50001, \$0)
            #    } else if (sn == "chrM" || sn == "MT") {
            #    printf("%05d\\t\\t%s\\n", 50002, \$0)
            #    } else {
            #    # fallback: sort these lexicographically after numeric contigs
            #    printf("100000\\t%s\\t%s\\n", sn, \$0)
            #    }
            #}
            #' "contigs.txt" > "contigs_keys.txt"

            # build new header: sort ##CONTIG lines
            #sort -t\$'\t' -k1,1n -k2,2 "contigs_keys.txt" | cut -f3- > "contigs_sorted.txt"
            #cat "start.txt" > "new_header.txt"
            #cat "contigs_sorted.txt" >> "new_header.txt"
            #cat "end.txt" >> "new_header.txt"

            # write sorted BAM file
            #rehead_name="\${base_name}_rehead.bam"
            #sort_name="\${base_name}_sorted.bam"
            #samtools reheader new_header.txt \$file > \$rehead_name
            #samtools sort \$rehead_name > \$sort_name

            ${params.conda}/bamixchecker/bin/picard AddOrReplaceReadGroups \\
                I=\$file \\
                O=\${base_name}_rg.bam \\
                SORT_ORDER=coordinate \\
                RGID=1 \\
                RGLB=lib1 \\
                RGPL=illumina \\
                RGPU=unit1 \\
                RGSM=\$base_name \\
                CREATE_INDEX=True \\
                REFERENCE_SEQUENCE=${params.refGenome}/fasta/genome.fa

            echo "\${base_name}_rg.bam" >> bam_list.txt
        done

        # Run BAMixChecker
        python ${params.BAMIXCHECKER}/BAMixChecker.py \\
            -l bam_list.txt \\
            -r ${params.refGenome}/fasta/genome.fa \\
            -o \${output_location} \\
            --OFFFileNameMatching \\
            -p 4
        """
}
