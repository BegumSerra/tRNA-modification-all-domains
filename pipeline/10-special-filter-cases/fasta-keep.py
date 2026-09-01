#!/usr/bin/env python3
"""
Keep sequences from a FASTA file that match IDs in a list file.
"""

import sys


def load_ids(id_file):
    """Load IDs from a text file (one ID per line)."""
    ids = set()
    with open(id_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


def filter_fasta_keep(input_fasta, id_file, output_fasta):
    """
    Keep only sequences whose IDs are in the ID list.
    
    Args:
        input_fasta: Path to input FASTA file
        id_file: Path to file with IDs to keep (one per line)
        output_fasta: Path to output FASTA file
    """
    # Load IDs to keep
    ids_to_keep = load_ids(id_file)
    print(f"Loaded {len(ids_to_keep)} IDs to keep from {id_file}")
    
    kept_count = 0
    total_count = 0
    
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        keep_sequence = False
        current_header = None
        
        for line in infile:
            if line.startswith('>'):
                total_count += 1
                # Extract the ID (first word after >)
                seq_id = line[1:].split()[0]
                
                # Check if this ID should be kept
                if seq_id in ids_to_keep:
                    keep_sequence = True
                    kept_count += 1
                    outfile.write(line)
                else:
                    keep_sequence = False
            else:
                if keep_sequence:
                    outfile.write(line)
    
    print(f"Total sequences in input: {total_count}")
    print(f"Sequences kept: {kept_count}")
    print(f"Sequences removed: {total_count - kept_count}")
    print(f"Output written to: {output_fasta}")
    
    # Check for IDs in list but not found in FASTA
    # This requires a second pass or tracking during first pass
    return kept_count


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python fasta-keep.py <input_fasta> <id_list_file> <output_fasta>")
        print("  Keep only sequences with IDs in the list file")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    id_file = sys.argv[2]
    output_fasta = sys.argv[3]
    
    filter_fasta_keep(input_fasta, id_file, output_fasta)
