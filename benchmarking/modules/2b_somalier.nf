#!/usr/bin/env nextflow

process SOMALIER_EXTRACT {
    conda "${params.conda}/somalier"
    publishDir "${params.outdir}/2b_somalier/extract", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2

    input:
        path(bam)
        path(bam_index)

    output:
        path("extracted/*.somalier")

    script:
        """
        set -euo pipefail

        if [ -z "${params.snp_vcf}" ]
        then
            echo "ERROR: params.snp_vcf is empty. Please set --snp_vcf /path/to/sites.hg38.vcf.gz" >&2
            exit 1
        fi

        mkdir -p extracted
        mkdir -p extracted/tmp
        # Keep temporary files task-local to avoid shared /scratch temp contention.
        export TMPDIR="\$PWD/extracted/tmp"

        sites_vcf="${params.snp_vcf}"

        # Somalier requires matching contig names between sites VCF and FASTA.
        fasta_contig="\$(grep -m1 '^>' "${params.refGenome}/fasta/genome.fa" | sed 's/^>//' | awk '{print \$1}')"
        # Avoid head/awk early-exit with pipefail (can surface as exit 141/SIGPIPE on HPC).
        vcf_contig="\$(bcftools query -f '%CHROM\n' "${params.snp_vcf}" | sed -n '1p')"

        if [[ -z "\$fasta_contig" || -z "\$vcf_contig" ]]; then
            echo "ERROR: Unable to determine contig naming style from FASTA or sites VCF." >&2
            exit 1
        fi

        if [[ "\$fasta_contig" == chr* && "\$vcf_contig" != chr* ]]; then
            echo "Detected contig mismatch: FASTA uses chr-prefix but sites VCF does not. Harmonizing sites VCF..."
            {
                for i in \$(seq 1 22); do
                    printf '%s\tchr%s\n' "\$i" "\$i"
                done
                printf 'X\tchrX\nY\tchrY\nMT\tchrM\nM\tchrM\n'
            } > extracted/contig_map.tsv
            bcftools annotate --rename-chrs extracted/contig_map.tsv -Oz -o extracted/sites.harmonized.unsorted.vcf.gz "${params.snp_vcf}"
            bcftools sort -Oz -o extracted/sites.harmonized.vcf.gz extracted/sites.harmonized.unsorted.vcf.gz
            tabix -f -p vcf extracted/sites.harmonized.vcf.gz
            sites_vcf="extracted/sites.harmonized.vcf.gz"
        elif [[ "\$fasta_contig" != chr* && "\$vcf_contig" == chr* ]]; then
            echo "Detected contig mismatch: FASTA does not use chr-prefix but sites VCF does. Harmonizing sites VCF..."
            {
                for i in \$(seq 1 22); do
                    printf 'chr%s\t%s\n' "\$i" "\$i"
                done
                printf 'chrX\tX\nchrY\tY\nchrM\tMT\nchrMT\tMT\n'
            } > extracted/contig_map.tsv
            bcftools annotate --rename-chrs extracted/contig_map.tsv -Oz -o extracted/sites.harmonized.unsorted.vcf.gz "${params.snp_vcf}"
            bcftools sort -Oz -o extracted/sites.harmonized.vcf.gz extracted/sites.harmonized.unsorted.vcf.gz
            tabix -f -p vcf extracted/sites.harmonized.vcf.gz
            sites_vcf="extracted/sites.harmonized.vcf.gz"
        fi

        # Ensure sites are compatible with the selected FASTA.
        # This drops records whose REF does not match the FASTA (e.g., FASTA base N).
        bcftools norm \
            --check-ref x \
            --fasta-ref "${params.refGenome}/fasta/genome.fa" \
            -Oz \
            -o extracted/sites.refchecked.unsorted.vcf.gz \
            "\$sites_vcf"
        bcftools sort -Oz -o extracted/sites.refchecked.vcf.gz extracted/sites.refchecked.unsorted.vcf.gz
        tabix -f -p vcf extracted/sites.refchecked.vcf.gz

        # Somalier only supports single-nucleotide variants for direct genotyping from BAM.
        bcftools view \
            -m2 -M2 \
            -v snps \
            -Oz \
            -o extracted/sites.somalier.vcf.gz \
            extracted/sites.refchecked.vcf.gz
        tabix -f -p vcf extracted/sites.somalier.vcf.gz
        sites_vcf="extracted/sites.somalier.vcf.gz"

        n_sites="\$(bcftools view -H "\$sites_vcf" | wc -l)"
        if [[ "\$n_sites" -eq 0 ]]; then
            echo "ERROR: No compatible biallelic SNP sites remain after harmonization/ref-check filtering." >&2
            echo "Check that params.snp_vcf and params.refGenome come from the same reference build." >&2
            exit 1
        fi

        # Somalier warns/limits to 65535 sites. Subset deterministically to reduce memory and improve stability.
        max_sites=65000
        if [[ "\$n_sites" -gt "\$max_sites" ]]; then
            step="\$(( (n_sites + max_sites - 1) / max_sites ))"
            bcftools view -h "\$sites_vcf" > extracted/sites.subset.header.vcf
            bcftools view -H "\$sites_vcf" | awk -v step="\$step" 'NR % step == 1' > extracted/sites.subset.body.vcf
            cat extracted/sites.subset.header.vcf extracted/sites.subset.body.vcf | bgzip -c > extracted/sites.somalier.subset.vcf.gz
            tabix -f -p vcf extracted/sites.somalier.subset.vcf.gz
            sites_vcf="extracted/sites.somalier.subset.vcf.gz"
            n_sites="\$(bcftools view -H "\$sites_vcf" | wc -l)"
            echo "Downsampled Somalier sites to \$n_sites loci (from >\$max_sites)."
        fi

        # Somalier requires either an @RG header in the BAM or SOMALIER_SAMPLE_NAME.
        # Derive a stable sample name from the BAM file name for BAMs without read-group tags.
        sample_name="\$(basename "${bam}" .bam)"
        export SOMALIER_SAMPLE_NAME="\$sample_name"

        somalier extract \
            --sites "\$sites_vcf" \
            --fasta "${params.refGenome}/fasta/genome.fa" \
            -d extracted \
            ${bam}
        """
}

process SOMALIER_RELATE {
    conda "${params.conda}/somalier"
    publishDir "${params.outdir}/2b_somalier", mode: 'copy'
    errorStrategy 'ignore'

    input:
        path(somalier_files)
        val(mod1)
        val(mod2)

    output:
        path("${params.dataset}/**/somalier/somalier*.tsv")
        path("${params.dataset}/**/somalier/somalier*.html"), optional: true

    script:
        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}"
        elif [ "${params.dataset}" == "hgsoc" ]
        then
            output_location="${params.dataset}/real_data/${mod1}_vs_${mod2}/ncells_null"
        else
            output_location="${params.dataset}/real_data/ncells_null"
        fi
        mkdir -p "\${output_location}/somalier"

        somalier relate \
            --infer \
            --output-prefix "\${output_location}/somalier/somalier" \
            ${somalier_files}
        """
}
