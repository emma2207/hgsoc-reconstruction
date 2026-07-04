process SNP_REFERENCE_INDEX {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/0_snp_ref_index", mode: 'copy'

    input:
    path vcf_file

    output:
    path "common_snps.vcf.gz"    , emit: vcf
    path "common_snps.vcf.gz.tbi", emit: tbi
    path "versions.yml"          , emit: versions

    script:
    """
    echo "Staging, chromosome re-naming (if needed), sorting, and indexing SNP reference: ${vcf_file}"
    
    echo "Input SNP VCF contigs:"
    bcftools view -H ${vcf_file} | cut -f1 | sort -u | head -30 || true
    
    ##########################################################################################################
    # Create chromosome rename map for GRCh38 / Cell Ranger style references
    #
    # Input common SNP resources often use:
    #   1, 2, 3, ..., 22, X, Y, MT
    #
    # Cell Ranger GRCh38 references usually use:
    #   chr1, chr2, chr3, ..., chr22, chrX, chrY, chrM
    #
    # If we do not rename the SNP VCF contigs, bcftools mpileup -T will find zero matching sites.
    ##########################################################################################################

    cat > chr_rename.txt << EOF
1\tchr1
2\tchr2
3\tchr3
4\tchr4
5\tchr5
6\tchr6
7\tchr7
8\tchr8
9\tchr9
10\tchr10
11\tchr11
12\tchr12
13\tchr13
14\tchr14
15\tchr15
16\tchr16
17\tchr17
18\tchr18
19\tchr19
20\tchr20
21\tchr21
22\tchr22
X\tchrX
Y\tchrY
MT\tchrM
M\tchrM
EOF

    ##########################################################################################################
    # Rename chromosomes, then sort and compress
    ##########################################################################################################

    echo "Renaming chromosomes where applicable..."

    bcftools annotate \\
        --rename-chrs chr_rename.txt \\
        ${vcf_file} \\
    | bcftools sort \\
        -O z \\
        -o common_snps.vcf.gz

    ##########################################################################################################
    # Index renamed/sorted SNP reference
    ##########################################################################################################

    bcftools index -t common_snps.vcf.gz
    
     ##########################################################################################################
    # Diagnostics
    ##########################################################################################################

    echo "Output SNP VCF contigs after renaming:"
    bcftools view -H common_snps.vcf.gz | cut -f1 | sort -u | head -30 || true

    SNP_COUNT=\$(bcftools view -H common_snps.vcf.gz | wc -l)
    echo "Output SNP count: \${SNP_COUNT}"

    if [ "\${SNP_COUNT}" -eq 0 ]; then
        echo "ERROR: SNP reference contains zero records after processing."
        exit 1
    fi

    echo "SNP reference indexing completed successfully:"
    ls -lh common_snps.vcf.gz common_snps.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
