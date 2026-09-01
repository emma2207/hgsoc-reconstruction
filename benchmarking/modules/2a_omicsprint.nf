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
        def omicsprintCallRate = (params.omicsprint_call_rate != null) ? params.omicsprint_call_rate : 0.80
        def omicsprintCoverageRate = (params.omicsprint_coverage_rate != null) ? params.omicsprint_coverage_rate : 0.25

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
                printf "\\n"
            }' > "\${matrix_body_file}"

        if [ "${params.omicsprint_max_snps}" -gt 0 ]
        then
            head -n "${params.omicsprint_max_snps}" "\${matrix_body_file}" > "\${matrix_body_file}.subset"
            mv "\${matrix_body_file}.subset" "\${matrix_body_file}"
        fi

        {
            printf "snp_id"
            awk '{printf "\t%s", \$1} END {printf "\\n"}' "\${sample_file}"
        } > "\${matrix_file}"
        cat "\${matrix_body_file}" >> "\${matrix_file}"

        {
            printf '%s\n' 'suppressPackageStartupMessages(library(omicsPrint))'
            printf '%s\n' ''
            printf '%s\n' 'args <- commandArgs(trailingOnly = TRUE)'
            printf '%s\n' 'genotype_file <- args[1]'
            printf '%s\n' 'out_prefix <- args[2]'
            printf '%s\n' 'call_rate <- ${omicsprintCallRate}'
            printf '%s\n' 'coverage_rate <- ${omicsprintCoverageRate}'
            printf '%s\n' ''
            printf '%s\n' 'geno_df <- read.delim(genotype_file, check.names = FALSE, stringsAsFactors = FALSE)'
            printf '%s\n' 'if (nrow(geno_df) == 0) {'
            printf '%s\n' '    stop("No SNP rows available for OMICSPRINT after filtering.")'
            printf '%s\n' '}'
            printf '%s\n' ''
            printf '%s\n' 'rownames(geno_df) <- geno_df[[1]]'
            printf '%s\n' 'x <- as.matrix(geno_df[, -1, drop = FALSE])'
            printf '%s\n' 'storage.mode(x) <- "numeric"'
            printf '%s\n' ''
            printf '%s\n' 'polymorphic <- apply(x, 1, function(row) {'
            printf '%s\n' '    length(unique(stats::na.omit(row))) > 1'
            printf '%s\n' '})'
            printf '%s\n' 'x <- x[polymorphic, , drop = FALSE]'
            printf '%s\n' ''
            printf '%s\n' 'if (nrow(x) < 2) {'
            printf '%s\n' '    stop("Need at least 2 polymorphic SNPs for OMICSPRINT.")'
            printf '%s\n' '}'
            printf '%s\n' 'if (ncol(x) < 2) {'
            printf '%s\n' '    stop("Need at least 2 samples for OMICSPRINT.")'
            printf '%s\n' '}'
            printf '%s\n' ''
            printf '%s\n' 'colnames(x) <- colnames(geno_df)[-1]'
            printf '%s\n' ''
            printf '%s\n' 'write_empty_outputs <- function(msg) {'
            printf '%s\n' '    plot_msg <- "No IBS scatter generated: insufficient high-quality samples after OMICSPRINT pruning."'
            printf '%s\n' '    message(msg)'
            printf '%s\n' '    write.table('
            printf '%s\n' '        data.frame(message = msg, stringsAsFactors = FALSE),'
            printf '%s\n' '        paste0(out_prefix, "_allele_sharing.tsv"),'
            printf '%s\n' '        sep = "\\t", quote = FALSE, row.names = FALSE'
            printf '%s\n' '    )'
            printf '%s\n' '    write.table('
            printf '%s\n' '        data.frame(message = msg, stringsAsFactors = FALSE),'
            printf '%s\n' '        paste0(out_prefix, "_mismatches.tsv"),'
            printf '%s\n' '        sep = "\\t", quote = FALSE, row.names = FALSE'
            printf '%s\n' '    )'
            printf '%s\n' '    grDevices::pdf(paste0(out_prefix, "_ibs_scatter.pdf"), width = 8, height = 6)'
            printf '%s\n' '    graphics::plot.new()'
            printf '%s\n' '    graphics::title(main = "omicsPrint IBS scatter unavailable")'
            printf '%s\n' '    graphics::text(0.5, 0.58, labels = plot_msg, cex = 0.9)'
            printf '%s\n' '    graphics::text(0.5, 0.42, labels = "See *_allele_sharing.tsv for the detailed reason.", cex = 0.8)'
            printf '%s\n' '    grDevices::dev.off()'
            printf '%s\n' '}'
            printf '%s\n' ''
            printf '%s\n' 'write_outputs <- function(relation_data) {'
            printf '%s\n' '    write.table('
            printf '%s\n' '        relation_data,'
            printf '%s\n' '        paste0(out_prefix, "_allele_sharing.tsv"),'
            printf '%s\n' '        sep = "\\t",'
            printf '%s\n' '        quote = FALSE,'
            printf '%s\n' '        row.names = FALSE'
            printf '%s\n' '    )'
            printf '%s\n' ''
            printf '%s\n' '    mismatches <- inferRelations(relation_data)'
            printf '%s\n' '    write.table('
            printf '%s\n' '        mismatches,'
            printf '%s\n' '        paste0(out_prefix, "_mismatches.tsv"),'
            printf '%s\n' '        sep = "\\t",'
            printf '%s\n' '        quote = FALSE,'
            printf '%s\n' '        row.names = FALSE'
            printf '%s\n' '    )'
            printf '%s\n' ''
            printf '%s\n' '    grDevices::pdf(paste0(out_prefix, "_ibs_scatter.pdf"), width = 8, height = 6)'
            printf '%s\n' '    graphics::plot(relation_data\$mean, relation_data\$var,'
            printf '%s\n' '                   xlab = "IBS mean", ylab = "IBS variance", pch = 16,'
            printf '%s\n' '                   main = "omicsPrint allele sharing")'
            printf '%s\n' '    grDevices::dev.off()'
            printf '%s\n' '}'
            printf '%s\n' ''
            printf '%s\n' 'attempts <- unique(rbind('
            printf '%s\n' '    data.frame(callRate = call_rate, coverageRate = coverage_rate),'
            printf '%s\n' '    data.frame(callRate = 0.60, coverageRate = 0.10),'
            printf '%s\n' '    data.frame(callRate = 0.50, coverageRate = 0.00),'
            printf '%s\n' '    data.frame(callRate = 0.30, coverageRate = 0.00)'
            printf '%s\n' '))'
            printf '%s\n' ''
            printf '%s\n' 'result <- FALSE'
            printf '%s\n' 'last_msg <- "OMICSPRINT failed for unknown reason."'
            printf '%s\n' 'for (i in seq_len(nrow(attempts))) {'
            printf '%s\n' '    cr <- attempts\$callRate[i]'
            printf '%s\n' '    cov <- attempts\$coverageRate[i]'
            printf '%s\n' '    message(paste0("Trying OMICSPRINT with callRate=", cr, ", coverageRate=", cov, " ..."))'
            printf '%s\n' '    ok <- tryCatch({'
            printf '%s\n' '        relation_data <- alleleSharing(x, verbose = TRUE, callRate = cr, coverageRate = cov)'
            printf '%s\n' '        write_outputs(relation_data)'
            printf '%s\n' '        TRUE'
            printf '%s\n' '    }, error = function(e) {'
            printf '%s\n' '        last_msg <<- paste0("OMICSPRINT failed after filtering/pruning (callRate=", cr, ", coverageRate=", cov, "): ", conditionMessage(e))'
            printf '%s\n' '        message(last_msg)'
            printf '%s\n' '        FALSE'
            printf '%s\n' '    })'
            printf '%s\n' '    if (ok) {'
            printf '%s\n' '        result <- TRUE'
            printf '%s\n' '        break'
            printf '%s\n' '    }'
            printf '%s\n' '}'
            printf '%s\n' ''
            printf '%s\n' 'if (!result) {'
            printf '%s\n' '    write_empty_outputs(last_msg)'
            printf '%s\n' '    message("Wrote placeholder OMICSPRINT outputs due to sparse post-pruning data.")'
            printf '%s\n' '}'
        } > run_omicsprint.R

        Rscript run_omicsprint.R "\${matrix_file}" "\${omicsprint_prefix}"
        """.stripIndent()
}