#!/bin/sh

# Check whether all files are present for each dataset and data type, and report any missing accessions.

datasets=("central_nervous_system_tumor" "central_nervous_system_tumor" "hgsoc_ariel" "high_grade_glioma" "high_grade_glioma" "low_grade_glioma" "wilms_tumor" "wilms_tumor")
data_types=("bulk" "single-nucleus" "all" "bulk" "single-cell" "all" "bulk" "single-nucleus")

echo "Checking for missing downloads..." > "missing_accessions.txt"

# Track down any failed downloads
for i in "${!datasets[@]}"; do
    dataset=${datasets[$i]}
    data_type=${data_types[$i]}
    fastq_location="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/${dataset}/${data_type}"
    echo "Checking dataset: $dataset, data type: $data_type"
    missing_accessions=()
    for accession in $(cat "../data/accession_lists/SRR_Acc_List_${data_type}_${dataset}.txt"); do
        if ! ls ${fastq_location}/${accession}*.fastq.gz 1> /dev/null 2>&1; then
            missing_accessions+=("$accession")
        fi
    done

    if [ ${#missing_accessions[@]} -eq 0 ]; then
        echo "All files for dataset: $dataset, data type: $data_type are present." >> "missing_accessions.txt"
    else
        echo "Missing files for dataset: $dataset, data type: $data_type:" >> "missing_accessions.txt"
        for acc in "${missing_accessions[@]}"; do
            echo "$acc" >> "missing_accessions.txt"
        done
    fi
done