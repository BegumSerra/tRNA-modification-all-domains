# Copyright (c) 2024, Begüm Serra Büyüktarakçı
# All rights reserved.

from Bio import SeqIO
import sys

input_file = sys.argv[1]
reference_fasta = sys.argv[2]
output_file = sys.argv[3]

def count_dash_line(sequence):
    return sequence.count('-')

def filter_and_match_sequences(input_file, reference_fasta, output_file):
    # Code 1: Filter sequences based on dash line length
    records = list(SeqIO.parse(input_file, "fasta"))
    filtered_records = [record for record in records if count_dash_line(record.seq) <= len(record) / 2]

    # Code 2: Match and filter sequences from the reference FASTA
    sequence_ids_to_keep = {record.id for record in filtered_records}
    reference_records = [record for record in SeqIO.parse(reference_fasta, "fasta") if record.id in sequence_ids_to_keep]

    # Write the matched and filtered sequences to a new FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(reference_records, output_handle, "fasta")

# Example usage:
filter_and_match_sequences(input_file, reference_fasta, output_file)
