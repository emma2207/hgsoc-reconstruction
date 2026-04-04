#!/usr/bin/env bash
set -euo pipefail

module load bcftools

in_dir="intersection_results"
out_dir="intersection_results_fixed"
mkdir -p "$out_dir"

shopt -s nullglob
files=( "$in_dir"/*.vcf.gz )

if [ ${#files[@]} -eq 0 ]; then
  echo "No .vcf.gz files found in $in_dir"
  exit 1
fi

for f in "${files[@]}"; do
  base=$(basename "$f")
  out_vcf="$out_dir/$base"
  hdr_tmp=$(mktemp)

  echo "Fixing $base"

  bcftools view -h "$f" 2>/dev/null | \
  sed -E 's/^##FORMAT=<ID=PL,Number=[^,]+,/##FORMAT=<ID=PL,Number=G,/' | \
  sed -E 's/^##INFO=<ID=DP4,Number=[^,]+,/##INFO=<ID=DP4,Number=4,/' \
  > "$hdr_tmp"

  bcftools reheader -h "$hdr_tmp" -o "$out_vcf" "$f"
  bcftools index -f "$out_vcf"

  rm -f "$hdr_tmp"

  echo "Header check for $base:"
  bcftools view -h "$out_vcf" | grep -E '^##FORMAT=<ID=PL|^##INFO=<ID=DP4'
  echo
done

echo "Done. Fixed files are in $out_dir"