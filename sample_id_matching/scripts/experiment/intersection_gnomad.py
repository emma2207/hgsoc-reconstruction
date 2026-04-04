import sys
import hail as hl

hl.init()

input_file = sys.argv[1]
filename = input_file.split("/")[-1]
sample_name = filename.split(".")[0]

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

mt.write(f"intersection_{sample_name}.mt", overwrite=True)
print("The end!")