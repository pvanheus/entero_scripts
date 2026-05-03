#!/usr/bin/env python3

import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat metadata for Clockor2")
    parser.add_argument(
        "metadata_filename",
        type=argparse.FileType(),
        help="Metadata file with Accession,Country,Collection_Date",
    )
    parser.add_argument(
        "output_file",
        type=argparse.FileType("w"),
        nargs="?",
        default=sys.stdout,
        help="Output file with extra columns",
    )
    args = parser.parse_args()

    line = next(args.metadata_filename)
    args.output_file.write("tip,date\n")
    for line in args.metadata_filename:
        accession, country, collection_date = line.strip().split(",")
        args.output_file.write(f"{accession},{collection_date}\n")
