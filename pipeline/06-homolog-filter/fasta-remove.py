#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
remove_file = sys.argv[2] # Input wanted file, one gene name per line
result_file = sys.argv[3] # Output fasta file

remove = set()
with open(remove_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            remove.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

with open(result_file, "w") as f:
    for seq in fasta_sequences:
        nam = seq.id
        nuc = str(seq.seq)
        if nam not in remove and len(nuc) > 0:
            SeqIO.write([seq], f, "fasta")

num_fasta = len([1 for line in open(fasta_file) if line.startswith(">")])
print("Number of sequences in the original file = ", num_fasta)

num_result = len([1 for line in open(result_file) if line.startswith(">")])
print("Number of sequences in the final file = ", num_result)

