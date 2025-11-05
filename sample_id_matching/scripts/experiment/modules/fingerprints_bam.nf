process CROSSCHECK_FINGERPRINTS_BAM {
    tag "${sample_id}"
    
    publishDir "${params.outdir}/fingerprints", mode: 'copy'
    
    input:
        tuple val(sample_id), path(bam)
        path haplotype_map
        val crosscheck_by // SAMPLE, LIBRARY, FILE
    
    output:
        tuple val(sample_id), path("${sample_id}.crosscheck_metrics.txt"), emit: metrics
    
    script:
        def memory = task.memory ? "-Xmx${task.memory.toGiga()}G" : '-Xmx3G'
        """
        set -euo pipefail
        
        # Create input list file
        echo "${bam}" > bam_list.txt
        
        # Run Picard CrosscheckFingerprints
        picard ${memory} CrosscheckFingerprints \\
            INPUT=bam_list.txt \\
            HAPLOTYPE_MAP=${haplotype_map} \\
            OUTPUT=${sample_id}.crosscheck_metrics.txt \\
            CROSSCHECK_BY=${crosscheck_by} \\
            NUM_THREADS=${task.cpus} \\
            VALIDATION_STRINGENCY=SILENT \\
            EXIT_CODE_WHEN_NO_VALID_CHECKS=1
        """
}