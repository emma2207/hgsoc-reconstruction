process NGSCHECKMATE_BAM {
    tag "${sample_id}"
    
    publishDir "${params.outdir}/ngscheckmate", mode: 'copy'
    
    input:
        tuple val(sample_id), path(bam)
        path ref_fasta
        path snp_bed
    
    output:
        tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf
        tuple val(sample_id), path("output_corr_matrix.txt"), path("output_all_results.txt"), emit: results
        path "r_script.r.Rout"
        path "plots/*"
    
    script:
        """
        set -euo pipefail
        
        # Create ncm.conf file
        echo "REF=\\"${ref_fasta}\\"" > ncm.conf
        echo "SAMTOOLS=\\"samtools\\"" >> ncm.conf
        echo "BCFTOOLS=\\"bcftools\\"" >> ncm.conf
        
        # Create bam list file
        echo "${bam}" > bam_list.txt
        
        # Create output directory
        mkdir -p plots
        
        # Run NGSCheckMate
        python ${params.NGS_CHECKMATE_HOME}/ncm.py \\
            -B \\
            -l bam_list.txt \\
            -bed ${snp_bed} \\
            -O . \\
            -N ${sample_id} \\
            -f
        
        # Rename outputs to standardized names
        [ -f "${sample_id}_corr_matrix.txt" ] && mv "${sample_id}_corr_matrix.txt" output_corr_matrix.txt
        [ -f "${sample_id}_all_results.txt" ] && mv "${sample_id}_all_results.txt" output_all_results.txt
        [ -f "r_script.r.pdf" ] && mv r_script.r.pdf plots/sample_clustering.pdf
        [ -f "${sample_id}.matched.result.txt" ] && mv "${sample_id}.matched.result.txt" plots/matched_samples.txt
        """
}