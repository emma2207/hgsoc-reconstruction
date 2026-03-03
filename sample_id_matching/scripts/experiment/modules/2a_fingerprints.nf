#!/usr/bin/env nextflow

process CROSSCHECK_FINGERPRINTS {
    conda "${params.conda}/fingerprints"
    publishDir "${params.outdir}/2a_fingerprints", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf)
    
    output:
        path("${params.dataset}/*/*/*/crosscheck_metrics.txt")

    script:
        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${params.read_depth}"
        else
            output_location="${params.dataset}/real_data/ncells_null/read_depth_${params.read_depth}" 
        fi
        mkdir -p \$output_location

        # Sort variants in vcf files numerically
        for vcf_file in ${vcf}
        do
            # extract header and isolate ##contig lines
            bcftools view -h "\$vcf_file" > "header.txt"
            grep -v '^##contig' "header.txt" | grep '^##' > "meta.txt"
            grep '^##contig' "header.txt" > "contigs.txt"
            grep '^#CHROM' "header.txt" > "colheader.txt"

            # produce sortable keys for contigs: numeric chrN -> N, chrX->50000, chrY->50001, chrM/MT->50002, others->100000+name
            awk '
            {
                if (match(\$0, /ID=[^,>]+/)) {
                sn = substr(\$0, RSTART + 3, RLENGTH - 3)
                } else {
                sn = ""
                }
                if (sn ~ /^chr[0-9]+\$/) {
                n = substr(sn,4) + 0
                printf("%05d\\t\\t%s\\n", n, \$0)
                } else if (sn == "chrX") {
                printf("%05d\\t\\t%s\\n", 50000, \$0)
                } else if (sn == "chrY") {
                printf("%05d\\t\\t%s\\n", 50001, \$0)
                } else if (sn == "chrM" || sn == "MT") {
                printf("%05d\\t\\t%s\\n", 50002, \$0)
                } else {
                # fallback: sort these lexicographically after numeric contigs
                printf("100000\\t%s\\t%s\\n", sn, \$0)
                }
            }
            ' "contigs.txt" > "contigs_keys.txt"

            # build new header: sort ##CONTIG lines
            sort -t\$'\t' -k1,1n -k2,2 "contigs_keys.txt" | cut -f3- > "contigs_sorted.txt"
            cat "meta.txt" > "new_header.txt"
            cat "contigs_sorted.txt" >> "new_header.txt"
            cat "colheader.txt" >> "new_header.txt"

            # write sorted VCF file
            base_name=\$(basename "\$vcf_file")
            base_name=\${base_name%.vcf.gz}
            base_name=\${base_name%.vcf}
            rehead_name="\${base_name}_rehead.vcf"
            sort_name="\${base_name}_sorted.vcf"
            bcftools reheader -h new_header.txt \$vcf_file > \$rehead_name
            bcftools sort \$rehead_name > \$sort_name

            if ! [[ \$vcf_file == *"2507"* && \$vcf_file == *"bulk_diss_polyA"* ]]
            then
                echo "\$sort_name" >> \${output_location}/vcf_list.txt
            fi
        done
        
        # Run Picard CrosscheckFingerprints
        ${params.conda}/fingerprints/bin/picard CrosscheckFingerprints \\
            INPUT=\${output_location}/vcf_list.txt \\
            HAPLOTYPE_MAP=${params.projectDir}/modules/hg38_chr.map \\
            OUTPUT=\${output_location}/crosscheck_metrics.txt \\
            CROSSCHECK_BY=FILE \\
            NUM_THREADS=4 \\
            VALIDATION_STRINGENCY=LENIENT \\
            EXIT_CODE_WHEN_MISMATCH=0    
        """
}