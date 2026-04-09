import sys
import hail as hl

hl.init()

input_file = sys.argv[1]
filename = input_file.split("/")[-1]
sample_name = filename.split(".")[0]
output_location = sys.argv[2]

mt = hl.import_vcf(
    input_file,
    force_bgz=True,
    reference_genome="GRCh38",
    skip_invalid_loci=True,
)
mt = hl.split_multi_hts(mt)

sites = hl.read_table(
    "gs://gcp-public-data--gnomad/release/4.1.1/ht/exomes/gnomad.exomes.v4.1.1.sites.ht"
)

# Keep only rows where (locus, alleles) exists in gnomAD sites
mt = mt.filter_rows(hl.is_defined(sites[mt.row_key]))

# Run variant QC to compute allele frequencies
mt = hl.variant_qc(mt)

# Filter rows based on MAF (e.g., > 1% / 0.01)
mt_filtered = mt.filter_rows(mt.variant_qc.AF[0] > 0.01)

hl.export_vcf(mt_filtered, f"{output_location}/{sample_name}.vcf.bgz")

print("The end!")
