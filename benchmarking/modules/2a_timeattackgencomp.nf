#!/usr/bin/env nextflow

include { buildReadDepthContext } from './helpers/read_depth_utils'

process TIMEATTACKGENCOMP {
    conda "${params.conda}/TimeAttackGenComp"
    publishDir "${params.outdir}/2a_timeattackgencomp", mode: 'copy'
    errorStrategy 'ignore'

    input:
        path(vcf)
        path(vcf_index)
        val(modalities)

    output:
        path("${params.dataset}/**/timeattackgencomp/timeattackgencomp.snv.out.txt")
        path("${params.dataset}/**/timeattackgencomp/timeattackgencomp.raw.out")
        path("${params.dataset}/**/timeattackgencomp/timeattackgencomp.pdf"), optional: true

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

        if [ "${params.pseudobulk}" == "true" ]
        then
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${pseudobulkReadDepth}"
        else
            output_location="${params.dataset}/${realDataNcellsPath}/read_depth_${modalityReadDepthTag}"
        fi
        mkdir -p "\${output_location}"

        if [ -z "${params.TIMEATTACKGENCOMP}" ]
        then
            echo "ERROR: params.TIMEATTACKGENCOMP is empty. Please set it to the TimeAttackGenComp directory." >&2
            exit 1
        fi

        if [ -z "${params.timeattackgencomp_target_bed}" ]
        then
            echo "ERROR: params.timeattackgencomp_target_bed is empty. Please set --timeattackgencomp_target_bed." >&2
            exit 1
        fi

        target_bed_path="${params.timeattackgencomp_target_bed}"
        if [[ "\$target_bed_path" == *.gz ]]
        then
            target_bed_path="\${output_location}/targets_for_timeattackgencomp.bed"
            gzip -dc "${params.timeattackgencomp_target_bed}" > "\$target_bed_path"
        fi

        snv_files=()

        for file in ${vcf}
        do
            base_name=\$(basename "\$file")
            base_name=\${base_name%.vcf.bgz}
            base_name=\${base_name%.vcf.gz}
            base_name=\${base_name%.vcf}

            sample_snv="\${output_location}/\${base_name}.smp.snv"
            bcftools query -R "\$target_bed_path" -f '%CHROM\t%POS\t[%TGT]\n' "\$file" \
                | awk 'BEGIN { OFS = "\t" } { if (\$3 != "./.") { split(\$3, alleles, "/"); if (alleles[1] > alleles[2]) { \$3 = alleles[2] alleles[1] } else { \$3 = alleles[1] alleles[2] } } print }' \
                > "\$sample_snv"
            snv_files+=("\$sample_snv")
        done

        if [ "\${#snv_files[@]}" -eq 0 ]
        then
            echo "ERROR: No .smp.snv files were generated in \${output_location}." >&2
            exit 1
        fi

        perl "${params.TIMEATTACKGENCOMP}/compare_simple.pl" --raw_output "\${output_location}/timeattackgencomp.raw.out" --bed "\$target_bed_path" --missing './.' "\${snv_files[@]}" > "\${output_location}/timeattackgencomp.snv.out.txt"

        # compare_simple.pl computes all-vs-all but only prints the upper triangle of the
        # matrix, leaving the lower triangle blank. Mirror the upper triangle into the lower
        # triangle so every pairwise comparison is present regardless of input file order.
        snv_matrix="\${output_location}/timeattackgencomp.snv.out.txt"
        awk -F'\t' '
            NR == 1 {
                for (i = 2; i <= NF; i++) { header[i] = \$i }
                ncols = NF
                next
            }
            NF == ncols {
                nrows++
                row[nrows] = \$1
                for (i = 2; i <= NF; i++) { val[nrows, i] = \$i }
                next
            }
            { trailing[++ntrail] = \$0 }
            END {
                for (i = 1; i <= nrows; i++) {
                    for (j = 1; j <= nrows; j++) {
                        if (val[i, j + 1] == "" && val[j, i + 1] != "") {
                            val[i, j + 1] = val[j, i + 1]
                        }
                    }
                }
                printf "\t"
                for (i = 2; i <= ncols; i++) { printf "%s%s", header[i], (i < ncols ? "\t" : "\n") }
                for (i = 1; i <= nrows; i++) {
                    printf "%s", row[i]
                    for (j = 2; j <= ncols; j++) { printf "\t%s", val[i, j] }
                    printf "\n"
                }
                for (t = 1; t <= ntrail; t++) { print trailing[t] }
            }
        ' "\$snv_matrix" > "\${snv_matrix}.symmetrized" && mv "\${snv_matrix}.symmetrized" "\$snv_matrix"

        if [ -f "${params.TIMEATTACKGENCOMP}/heatmap.R" ]
        then
            R --vanilla < "${params.TIMEATTACKGENCOMP}/heatmap.R" --args "\${output_location}/timeattackgencomp.snv.out.txt" "\${output_location}/timeattackgencomp"
        fi
        """.stripIndent()
}
