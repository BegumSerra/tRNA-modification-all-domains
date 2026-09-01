#!/usr/bin/env python3
"""
Reorder protein columns in the presence/absence table according to functional order
"""

import sys
import pandas as pd

def main():
    if len(sys.argv) != 4:
        print("Usage: python reorder-proteins.py <input_tsv> <protein_order_file> <output_tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    order_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Read the input table
    df = pd.read_csv(input_file, sep='\t')
    
    # Define protein renaming
    protein_renames = {
        'MnmE': 'MnmE/GTPB3',
        'MnmG': 'MnmG/MTO1',
        'IscS': 'IscS/NFS1',
        'Tuc1': 'Tuc1/CTU1',
        'Elp3': 'ELP3'
    }
    
    # Proteins to remove
    proteins_to_remove = ['GTPB3', 'MTO1', 'NFS1', 'CTU1']
    
    # Get taxonomy columns (first 5 columns)
    taxonomy_cols = ['Domain', 'Phylum', 'Order', 'Family', 'Species']
    
    # Rename proteins
    rename_dict = {}
    for old_name, new_name in protein_renames.items():
        if old_name in df.columns:
            rename_dict[old_name] = new_name
            print(f"Renaming {old_name} to {new_name}")
    
    if rename_dict:
        df = df.rename(columns=rename_dict)
    
    # Remove specified proteins
    proteins_removed = []
    for protein in proteins_to_remove:
        if protein in df.columns:
            df = df.drop(columns=[protein])
            proteins_removed.append(protein)
    
    if proteins_removed:
        print(f"Removed proteins: {proteins_removed}")
    
    # Read the protein order
    with open(order_file, 'r') as f:
        protein_order_text = f.read().strip()
    
    # Split by tabs or spaces to get protein names
    protein_order = [p.strip() for p in protein_order_text.replace('\t', ' ').split() if p.strip()]
    
    # Get current protein columns (excluding taxonomy)
    current_protein_cols = [col for col in df.columns if col not in taxonomy_cols]
    
    # Reorder proteins based on the order file
    # Keep proteins that are in the order file in that order
    ordered_proteins = []
    for protein in protein_order:
        if protein in current_protein_cols:
            ordered_proteins.append(protein)
    
    # Add any remaining proteins that weren't in the order file
    for protein in current_protein_cols:
        if protein not in ordered_proteins:
            ordered_proteins.append(protein)
    
    # Create the final column order: taxonomy + ordered proteins
    final_columns = taxonomy_cols + ordered_proteins
    
    # Reorder the dataframe
    df_reordered = df[final_columns]
    
    # Save to output
    df_reordered.to_csv(output_file, sep='\t', index=False)
    
    print(f"Reordered {len(ordered_proteins)} protein columns")
    print(f"Output saved to: {output_file}")

if __name__ == '__main__':
    main()