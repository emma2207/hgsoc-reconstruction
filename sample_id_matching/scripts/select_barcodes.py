import argparse
import os
import pysam
import random

if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d", "--dataset", dest="dataset", help="Dataset name", required=True, type=str
    )
    parser.add_argument(
        "-s",
        "--sample",
        dest="sample",
        help="Sample name",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-b",
        "--bam",
        dest="bam_path",
        help="BAM file to read barcodes from",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-n",
        "--n_barcodes",
        dest="n_barcodes",
        type=int,
        default=1000,
        help="number of barcodes to select",
        required=True,
    )
    parser.add_argument(
        "-j",
        "--n_pseudobulks",
        dest="n_pseudobulks",
        type=int,
        default=1,
        help="number of pseudobulks to generate from a single file",
        required=False,
    )

    args = parser.parse_args()

    # Extract unique cell barcodes from BAM file
    bam_file = pysam.AlignmentFile(args.bam_path, "rb")

    cell_barcodes = [read.get_tag("CB") for read in bam_file if read.has_tag("CB")]
    cell_barcodes = [b for b in cell_barcodes if len(str(b)) == 16]
    unique_barcodes = list(set(cell_barcodes))

    for i in range(args.n_pseudobulks):
        # Randomly sample barcodes
        random_barcodes = random.sample(unique_barcodes, args.n_barcodes)

        subset_dir = f"subset-bam/{args.dataset}/{args.sample}"
        os.makedirs(subset_dir, exist_ok=True)

        # Write to file
        with open(
            f"{subset_dir}/selected_barcodes_{i+1}.txt",
            "w",
        ) as f:
            for barcode in random_barcodes:
                f.write(f"{barcode}\n")
