#!/usr/bin/env python3
"""
Add sequences to fasta file based on ID list from a01 folder.
Retrieves sequences from BLAST database files.
"""

import sys
from Bio import SeqIO
from pathlib import Path
import os

def read_id_list(id_file):
    """Read protein IDs from text file"""
    ids = set()
    with open(id_file) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.add(line)
    return ids

def get_sequences_from_database(protein_ids, db_dir):
    """
    Retrieve sequences for given protein IDs from database files.
    
    Args:
        protein_ids: Set of protein IDs to retrieve
        db_dir: Directory containing archaea, bacteria, eukarya subdirectories
    
    Returns:
        Dictionary mapping protein ID to SeqRecord
    """
    sequences = {}
    found_ids = set()
    
    # Check archaea, bacteria, and eukarya databases
    for domain in ['archaea', 'bacteria', 'eukarya']:
        domain_dir = Path(db_dir) / domain
        if not domain_dir.exists():
            print(f"Warning: {domain_dir} does not exist")
            continue
        
        # Look for fasta/faa files in the domain directory
        for fasta_file in list(domain_dir.glob('*.fasta')) + list(domain_dir.glob('*.faa')):
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # Check if this ID is in our list
                # Extract just the accession (e.g., WP_011019471.1 from >WP_011019471.1 description)
                record_id = record.id.split()[0].split('|')[-1]
                
                if record_id in protein_ids and record_id not in found_ids:
                    sequences[record_id] = record
                    found_ids.add(record_id)
                    print(f"Found {record_id} in {fasta_file.name}")
    
    # Report missing IDs
    missing = protein_ids - found_ids
    if missing:
        print(f"\nWarning: Could not find {len(missing)} sequences:")
        for mid in sorted(list(missing)[:10]):  # Show first 10
            print(f"  {mid}")
        if len(missing) > 10:
            print(f"  ... and {len(missing) - 10} more")
    
    print(f"\nFound {len(sequences)} out of {len(protein_ids)} requested sequences")
    return sequences

def add_sequences_to_fasta(input_fasta, id_list_file, output_fasta, db_dir):
    """
    Add sequences from ID list to existing fasta file.
    
    Args:
        input_fasta: Original fasta file
        id_list_file: File containing IDs to add
        output_fasta: Output fasta file with added sequences
        db_dir: Database directory with domain subdirectories
    """
    # Read IDs to add
    ids_to_add = read_id_list(id_list_file)
    print(f"Read {len(ids_to_add)} IDs from {id_list_file}")
    
    # Get existing IDs from input fasta
    existing_ids = set()
    existing_records = []
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = record.id.split('|')[-1] if '|' in record.id else record.id
        existing_ids.add(record_id)
        existing_records.append(record)
    
    print(f"Input fasta has {len(existing_records)} sequences")
    
    # Filter out IDs that already exist
    new_ids = ids_to_add - existing_ids
    already_present = ids_to_add & existing_ids
    
    if already_present:
        print(f"{len(already_present)} IDs already in fasta file (will skip)")
    
    if not new_ids:
        print("No new sequences to add!")
        # Just copy input to output
        with open(output_fasta, 'w') as out_f:
            SeqIO.write(existing_records, out_f, 'fasta')
        return
    
    print(f"Adding {len(new_ids)} new sequences")
    
    # Retrieve sequences from database
    new_sequences = get_sequences_from_database(new_ids, db_dir)
    
    # Write combined fasta
    with open(output_fasta, 'w') as out_f:
        # Write existing sequences
        SeqIO.write(existing_records, out_f, 'fasta')
        # Write new sequences
        SeqIO.write(list(new_sequences.values()), out_f, 'fasta')
    
    print(f"\nWrote {len(existing_records) + len(new_sequences)} total sequences to {output_fasta}")
    print(f"  Original: {len(existing_records)}")
    print(f"  Added: {len(new_sequences)}")

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python add-seqs.py <input_fasta> <id_list_file> <output_fasta> <database_dir>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    id_list_file = sys.argv[2]
    output_fasta = sys.argv[3]
    db_dir = sys.argv[4]
    
    add_sequences_to_fasta(input_fasta, id_list_file, output_fasta, db_dir)
