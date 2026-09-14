#!/usr/bin/env nextflow

process CONPAIR {
    conda "${params.conda}/conpair"
    publishDir "${params.outdir}/2b_conpair", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        tuple val(sample_1), path(bam_1), val(sample_2), path(bam_2), val(mod_1), val(mod_2)
    
    output:
        path("${params.dataset}/**/${sample_1}_vs_${sample_2}_concordance.txt")
    
    script:
        def conpairMarkerBed = params.conpair_marker_bed ?: "${params.CONPAIR}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
        def conpairMarkerTxt = params.conpair_marker_txt ?: "${params.CONPAIR}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
        def conpairMinCov = params.conpair_min_cov ?: 3
        def conpairMinMapQual = params.conpair_min_map_qual ?: 0
        def conpairMinBaseQual = params.conpair_min_base_qual ?: 10

        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}"
        elif [ "${params.dataset}" == "hgsoc" ]
        then
            output_location="${params.dataset}/real_data/${mod_1}_vs_${mod_2}/ncells_null"
        else
            output_location="${params.dataset}/real_data/ncells_null"
        fi

        export CONPAIR_DIR=/projects/\${USER}/software/conpair
        export GATK_JAR=/projects/\${USER}/software/anaconda/envs/conpair/opt/gatk-3.8/GenomeAnalysisTK.jar
        export PYTHONPATH=/projects/\${USER}/software/conpair/modules

        conpair_marker_bed="${conpairMarkerBed}"
        conpair_marker_txt="${conpairMarkerTxt}"
        conpair_min_cov=${conpairMinCov}
        conpair_min_map_qual=${conpairMinMapQual}
        conpair_min_base_qual=${conpairMinBaseQual}
        conpair_normalize_pileup_py="${params.projectDir}/modules/helpers/conpair_normalize_pileup.py"
        pileup_1="${sample_1}_pileup.txt"
        pileup_1_raw="${sample_1}_pileup.raw.txt"
        pileup_2="${sample_2}_pileup.txt"
        pileup_2_raw="${sample_2}_pileup.raw.txt"
        output_file="\${output_location}/${sample_1}_vs_${sample_2}_concordance.txt"

        if [ ! -s "\$conpair_marker_bed" ]
        then
            echo "ERROR: Conpair BED marker file not found: \$conpair_marker_bed" >&2
            exit 1
        fi

        if [ ! -s "\$conpair_marker_txt" ]
        then
            echo "ERROR: Conpair TXT marker file not found: \$conpair_marker_txt" >&2
            exit 1
        fi

        if [ ! -s "\$conpair_normalize_pileup_py" ]
        then
            echo "ERROR: Conpair pileup normalizer not found: \$conpair_normalize_pileup_py" >&2
            exit 1
        fi

        normalize_and_sanitize_pileup() {
            local raw_pileup="\$1"
            local out_pileup="\$2"
            local tmp_pileup="\${out_pileup}.tmp"

            python "\$conpair_normalize_pileup_py" "\$raw_pileup" "\$tmp_pileup"

            # Keep only lines with aligned base and quality strings for Conpair parsing.
            python - "\$tmp_pileup" "\$out_pileup" <<'PY'
import sys

in_path, out_path = sys.argv[1], sys.argv[2]

with open(in_path, "r", encoding="utf-8", errors="replace") as fin, open(out_path, "w", encoding="utf-8") as fout:
    for line in fin:
        if not line.strip() or line.startswith("[REDUCE RESULT]"):
            continue

        fields = line.rstrip("\\n").split("\\t")
        if len(fields) < 5:
            fields = line.rstrip("\\n").split()
        if len(fields) < 5:
            continue

        chrom, pos, ref = fields[0], fields[1], fields[2].upper()
        bases = fields[3].upper()
        quals = fields[4]

        # Zip truncates to the shorter sequence to guarantee index-safe alignment.
        filtered_pairs = [(b, q) for b, q in zip(bases, quals) if b in {"A", "C", "G", "T"}]
        if not filtered_pairs:
            continue

        out_bases = "".join(b for b, _ in filtered_pairs)
        out_quals = "".join(q for _, q in filtered_pairs)

        if len(out_bases) != len(out_quals):
            continue

        fout.write(f"{chrom}\\t{pos}\\t{ref}\\t{out_bases}\\t{out_quals}\\n")
PY

            rm -f "\$tmp_pileup"

            if [ ! -s "\$out_pileup" ]
            then
                echo "ERROR: Sanitized Conpair pileup is empty: \$out_pileup" >&2
                return 1
            fi
        }

        fasta_contig="\$(grep -m1 '^>' "${params.refGenome}/fasta/genome.fa" | sed 's/^>//' | awk '{print \$1}')"
        marker_contig="\$(awk 'NF {print \$1; exit}' "\$conpair_marker_txt")"
        remove_chr_opt=""

        if [[ "\$fasta_contig" == chr* && "\$marker_contig" != chr* ]]
        then
            remove_chr_opt="--remove_chr_prefix"
        elif [[ "\$fasta_contig" != chr* && "\$marker_contig" == chr* ]]
        then
            echo "ERROR: Reference FASTA contigs do not use 'chr' prefix, but selected Conpair markers do: \$conpair_marker_txt" >&2
            echo "Set --conpair_marker_bed/--conpair_marker_txt to marker files that match your reference naming." >&2
            exit 1
        fi

        # Run CONPAIR for one explicit pair per task.
        python ${params.CONPAIR}/scripts/run_gatk_pileup_for_sample.py \
                -B ${bam_1} \
                -R ${params.refGenome}/fasta/genome.fa \
                -M \$conpair_marker_bed \
                \$remove_chr_opt \
                -O \${pileup_1_raw}

        normalize_and_sanitize_pileup "\${pileup_1_raw}" "\${pileup_1}"

        python ${params.CONPAIR}/scripts/run_gatk_pileup_for_sample.py \
                -B ${bam_2} \
                -R ${params.refGenome}/fasta/genome.fa \
                -M \$conpair_marker_bed \
                \$remove_chr_opt \
                -O \${pileup_2_raw}

        normalize_and_sanitize_pileup "\${pileup_2_raw}" "\${pileup_2}"

        pileup_contig_1="\$(awk 'NF {print \$1; exit}' \${pileup_1})"
        pileup_contig_2="\$(awk 'NF {print \$1; exit}' \${pileup_2})"
        if [[ -z "\$pileup_contig_1" || -z "\$pileup_contig_2" ]]
        then
            echo "ERROR: Empty pileup detected for one or both samples: \${pileup_1}, \${pileup_2}" >&2
            exit 1
        fi

        mkdir -p "\$(dirname "\${output_file}")"

        python ${params.CONPAIR}/scripts/verify_concordance.py \
            -T \${pileup_1} \
            -N \${pileup_2} \
            -M \$conpair_marker_txt \
            -C \$conpair_min_cov \
            -Q \$conpair_min_map_qual \
            -B \$conpair_min_base_qual \
            -O \${output_file}

        # Conpair exits successfully even when no markers pass filters and may skip writing output.
        # Retry once with permissive thresholds to recover concordance for low-depth datasets.
        if [ ! -s "\${output_file}" ]
        then
            echo "WARN: No concordance output for \${sample_1} vs \${sample_2} at C=\$conpair_min_cov,Q=\$conpair_min_map_qual,B=\$conpair_min_base_qual; retrying with C=1,Q=0,B=0" >&2
            python ${params.CONPAIR}/scripts/verify_concordance.py \
                -T \${pileup_1} \
                -N \${pileup_2} \
                -M \$conpair_marker_txt \
                -C 1 \
                -Q 0 \
                -B 0 \
                -O \${output_file}
        fi

        if [ ! -s "\${output_file}" ]
        then
            echo "ERROR: Conpair did not produce concordance output for \${sample_1} vs \${sample_2} even after fallback thresholds." >&2
            echo "Check pileup depth and marker/reference compatibility." >&2
            exit 1
        fi
        """
}