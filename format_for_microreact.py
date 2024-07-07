#!/usr/bin/env python3

import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add year,month,day columns for Microreact')
    parser.add_argument('metadata_filename', type=argparse.FileType(), help='Metadata file with Accession,Country,Collection_Date')
    parser.add_argument('output_file', type=argparse.FileType('w'), nargs='?', default=sys.stdout, help='Output file with extra columns')
    args = parser.parse_args()

    line = next(args.metadata_filename)
    line = line.strip() + ',year,month,day\n'
    args.output_file.write(line)
    for line in args.metadata_filename:
        accession, country, collection_date = line.strip().split(',')
        year = month = day = ''
        date_parts = collection_date.split('-')
        year = date_parts[0]
        if len(date_parts) > 1:
            month = date_parts[1]
        if len(date_parts) > 2:
            day = date_parts[2]
        args.output_file.write(f'{accession},{country},{collection_date},{year},{month},{day}\n')
