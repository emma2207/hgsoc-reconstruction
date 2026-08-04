#!/usr/bin/env nextflow

include { buildReadDepthContext } from './helpers/read_depth_utils'

process OMICSPRINT {
    conda "${params.conda}/omicsPrint"
    publishDir "${params.outdir}/2a_omicsprint", mode: 'copy'
    errorStrategy 'ignore'

    input:
        path(vcf)
        val(modalities)

    output:
        path("${params.dataset}/**/omicsprint/omicsprint_genotype_matrix.tsv")
        path("${params.dataset}/**/omicsprint/omicsprint_input_samples.txt")
        path("${params.dataset}/**/omicsprint/omicsprint_allele_sharing.tsv")
        path("${params.dataset}/**/omicsprint/omicsprint_mismatches.tsv")
        path("${params.dataset}/**/omicsprint/omicsprint_ibs_scatter.pdf")

    script:
        def readDepthContext = buildReadDepthContext(
            params.read_depth,
            params.dataset,
            params.pseudobulk,
            modalities
        )

        def modalityReadDepthTag = readDepthContext.modalityReadDepthTag
        def pseudobulkReadDepth = readDepthContext.pseudobulkReadDepth
        def realDataNcellsPath = readDepthContext.realDataNcellsPath

        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${pseudobulkReadDepth}"
        else
            output_location="${params.dataset}/${realDataNcellsPath}/read_depth_${modalityReadDepthTag}"
        fi
        mkdir -p "\${output_location}/omicsprint"

        sample_file="\${output_location}/omicsprint/omicsprint_input_samples.txt"
        matrix_file="\${output_location}/omicsprint/omicsprint_genotype_matrix.tsv"
        matrix_body_file="\${output_location}/omicsprint/omicsprint_genotype_matrix.body.tsv"
        omicsprint_prefix="\${output_location}/omicsprint/omicsprint"

        bcftools query -l "${vcf}" > "\${sample_file}"

        n_samples=\$(wc -l < "\${sample_file}")
        if [ "\${n_samples}" -lt 2 ]
        then
            echo "ERROR: OMICSPRINT requires at least 2 samples. Found \${n_samples}." >&2
            exit 1
        fi

        # Keep only bi-allelic SNPs and convert GT strings to numeric classes expected by omicsPrint.
        bcftools view -m2 -M2 -v snps "${vcf}" \
            | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
            | awk 'BEGIN { OFS = "\t" } {
                snp = \$1 ":" \$2 ":" \$3 ":" \$4
                printf "%s", snp
                for (i = 5; i <= NF; i++) {
                    gt = \$i
                    val = "NA"
                    if (gt == "0/0" || gt == "0|0") {
                        val = 1
                    } else if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
                        val = 2
                    } else if (gt == "1/1" || gt == "1|1") {
                        val = 3
                    }
                    printf OFS val
                }
                printf "\n"
            }' > "\${matrix_body_file}"

        if [ "${params.omicsprint_max_snps}" -gt 0 ]
        then
            head -n "${params.omicsprint_max_snps}" "\${matrix_body_file}" > "\${matrix_body_file}.subset"
            mv "\${matrix_body_file}.subset" "\${matrix_body_file}"
        fi

        {
            printf "snp_id"
            awk '{printf "\t%s", $1} END {printf "\n"}' "\${sample_file}"
        } > "\${matrix_file}"
        cat "\${matrix_body_file}" >> "\${matrix_file}"

        cat > run_omicsprint.R << 'EOF'
        suppressPackageStartupMessages(library(omicsPrint))

        args <- commandArgs(trailingOnly = TRUE)
        genotype_file <- args[1]
        out_prefix <- args[2]

        geno_df <- read.delim(genotype_file, check.names = FALSE, stringsAsFactors = FALSE)
        if (nrow(geno_df) == 0) {
            stop("No SNP rows available for OMICSPRINT after filtering.")
        }

        rownames(geno_df) <- geno_df[[1]]
        x <- as.matrix(geno_df[, -1, drop = FALSE])
        storage.mode(x) <- "numeric"

        polymorphic <- apply(x, 1, function(row) {
            length(unique(stats::na.omit(row))) > 1
        })
        x <- x[polymorphic, , drop = FALSE]

        if (nrow(x) < 2) {
            stop("Need at least 2 polymorphic SNPs for OMICSPRINT.")
        }
        if (ncol(x) < 2) {
            stop("Need at least 2 samples for OMICSPRINT.")
        }

        relation_data <- alleleSharing(x, verbose = TRUE)
        write.table(
            relation_data,
            paste0(out_prefix, "_allele_sharing.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        mismatches <- inferRelations(relation_data)
        write.table(
            mismatches,
            paste0(out_prefix, "_mismatches.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        grDevices::pdf(paste0(out_prefix, "_ibs_scatter.pdf"), width = 8, height = 6)
        graphics::plot(relation_data\$mean, relation_data\$var,
                       xlab = "IBS mean", ylab = "IBS variance", pch = 16,
                       main = "omicsPrint allele sharing")
        grDevices::dev.off()
        EOF

        Rscript run_omicsprint.R "\${matrix_file}" "\${omicsprint_prefix}"
        """
}