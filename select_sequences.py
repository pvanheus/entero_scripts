#!/usr/bin/env python3

import argparse
import re

from csv import DictReader,DictWriter
from typing import TextIO

from Bio import SeqIO

def select_vp1_segments(vadr_output: TextIO):
    for line in vadr_output:
        if line.startswith('#'):
            continue
        fields = line.strip().split()
        # the VP1 protein in the model is encoded by a 915 bp segment
        if fields[3] == 'PASS' and fields[6] == 'vp1' and int(fields[14]) >= 900:
            name = fields[1]
            start = int(fields[10])
            end = int(fields[11])
            yield name, start, end

def select_good_samples(metadata_file: TextIO) -> set[str]:
    date_re = re.compile(r'^\d{4}(-\d{2}){0,2}$')
    good_samples = {}
    for entry in DictReader(metadata_file):
        if entry['Country'] and date_re.match(entry['Collection_Date']):
            year = int(entry['Collection_Date'].split('-')[0])
            if year > 1988:
                good_samples[entry['Accession']] = entry
    return good_samples

def select_sequences(vadr_output: TextIO, metadata_file: TextIO, fasta: TextIO, metadata_output_file: TextIO, limit_to: TextIO):
    vp1_coordinates = {name: (start, end) for name, start, end in select_vp1_segments(vadr_output)}
    metadata_writer = None
    if metadata_file:
        good_samples = select_good_samples(metadata_file)
        if metadata_output_file:
            if good_samples:
                record_id = list(good_samples.keys())[0]
                field_names = list(good_samples[record_id].keys())
                metadata_writer = DictWriter(metadata_output_file, fieldnames=field_names)
                metadata_writer.writeheader()
    else:
        good_samples = None
    if limit_to:
        acceptable_names = set([line.strip() for line in limit_to])

    for record in SeqIO.parse(fasta, 'fasta'):
        if record.id in vp1_coordinates and (good_samples is None or record.id in good_samples) and (limit_to is None or record.id in acceptable_names):
            if metadata_writer:
                metadata_writer.writerow(good_samples[record.id])
            start, end = vp1_coordinates[record.id]
            record.description += f' vp1|{start}|{end}|{end-start+1}bp'
            yield record[start-1:end]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select sequences for further analysis')
    parser.add_argument('--metadata_output_file', type=argparse.FileType('w'), help='Output file for metadata from selected samples')
    parser.add_argument('--limit_to', type=argparse.FileType(), help='List of sequence names that restrict what should be selected')
    parser.add_argument('vadr_output', type=argparse.FileType(), help='.sgm file output from VADR')
    parser.add_argument('fasta', type=argparse.FileType(), help='fasta file containing sequences to be selected')
    parser.add_argument('metadata_file', type=argparse.FileType(), nargs='?', help='metadata file containing sample information (optional)')
    args = parser.parse_args()
    for record in select_sequences(args.vadr_output, args.metadata_file, args.fasta, args.metadata_output_file, args.limit_to):
        print(record.format('fasta'), end='')
