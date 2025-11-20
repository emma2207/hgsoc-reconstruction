#!/usr/bin/env nextflow

process CREATE_FINGERPRINT_MAP {
    tag "fingerprint map"
    conda "fingerprint_map"
    publishDir "${params.outdir}/fingerprints", mode: 'copy'
    
    input:
        path(vcf)
    
    output:
        path("hg38_chr.map"), emit: map
    
    script:
        """
        set -euo pipefail
        
        # Create fingerprint map
        echo "Fingerprint map"
        """
}