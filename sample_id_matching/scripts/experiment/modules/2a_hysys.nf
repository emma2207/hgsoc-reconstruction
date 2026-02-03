#!/usr/bin/env nextflow

process HYSYS {
    conda "${params.conda}/hysys"
    publishDir "${params.outdir}/2a_hysys", mode: 'copy'
    errorStrategy 'ignore' 
    
    input:
        path(vcf)
        val(modalities)
    
    output:
        path("${params.dataset}/*/concordance_output.txt")
        path("${params.dataset}/*/model_results.txt")
    
    script:
        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${params.read_depth}"
        else
            output_location="${params.dataset}/real_data/read_depth_${params.read_depth}" 
        fi
        mkdir -p \$output_location

        # Prep HYSYS input
        for file in ${vcf}
        do
            # get filename
            base_name=\$(basename "\$file")
            base_name=\${base_name%.vcf.gz}
            base_name=\${base_name%.vcf}

            bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' \
                \${file} -o \${output_location}/\${base_name}.vcf.gz

            echo "Successfully annotated \${file}"

            # Add AF field to VCFs and convert VCFs to required format 
            # (SNP_ID CHR POS VAF)
            bcftools +fill-tags \${output_location}/\${base_name}.vcf.gz -- -t AF | \
                bcftools query -f '%ID\t%CHROM\t%POS\t%AF\n' | \
            awk '\$4!="."' | sort -k1,1 > "\${output_location}/\${base_name}.snps"

            echo "Created \${base_name}.snps"

            # Add file-path to list
            for mod in ${modalities.join(' ')}
            do  
                if [[ "\$base_name" == *"\$mod"* ]]
                then
                    if ! [[ \$file == *"2507"* && \$file == *"bulk_diss_polyA"* ]]
                    then
                        echo "\${output_location}/\${base_name}.snps" >> \\
                            \${output_location}/sample_list_\${mod}.txt
                        echo "Added \${base_name} to sample_list_\${mod}.txt"
                    fi
                fi
            done
            
        done
        
        # HYSYS
        # Run concordance calculation
        if [ ${params.pseudobulk} == true ]
        then 
            ${params.HYSYS}/HaveYouSwappedYourSamples.sh conc -s \\
                \${output_location}/sample_list_${modalities[0]}.txt \\
                \${output_location}/sample_list_${modalities[0]}.txt \\
                \${output_location}/concordance_output.txt
        else
            ${params.HYSYS}/HaveYouSwappedYourSamples.sh conc -s \\
                \${output_location}/sample_list_${modalities[0]}.txt \\
                \${output_location}/sample_list_${modalities[1]}.txt \\
                \${output_location}/concordance_output.txt
        fi
        
        echo "Ran concordance"
            
        # Run model analysis for detecting swaps
        ${params.HYSYS}/HaveYouSwappedYourSamples.sh model \\
            \${output_location}/concordance_output.txt \\
            > \${output_location}/model_results.txt

        echo "Ran modelling"
        """
}