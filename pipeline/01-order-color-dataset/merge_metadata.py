#!/usr/bin/env python3
"""
Merge archaea and bacteria GTDB filtered metadata files.
Combines rows from both tables into a single merged metadata file.
"""

import pandas as pd
import argparse
from pathlib import Path


def merge_metadata(bacteria_file, archaea_file, output_file):
    """
    Merge bacteria and archaea metadata files by concatenating rows.
    
    Parameters:
    -----------
    bacteria_file : str
        Path to bacteria filtered metadata TSV file
    archaea_file : str
        Path to archaea filtered metadata TSV file
    output_file : str
        Path to save merged metadata TSV file
    """
    print(f"Loading bacteria metadata from {bacteria_file}...")
    bacteria_df = pd.read_csv(bacteria_file, sep='\t')
    print(f"  Bacteria assemblies: {len(bacteria_df)}")
    
    print(f"\nLoading archaea metadata from {archaea_file}...")
    archaea_df = pd.read_csv(archaea_file, sep='\t')
    print(f"  Archaea assemblies: {len(archaea_df)}")
    
    # Verify that both dataframes have the same columns
    bacteria_cols = set(bacteria_df.columns)
    archaea_cols = set(archaea_df.columns)
    
    if bacteria_cols != archaea_cols:
        print("\nWARNING: Column mismatch detected!")
        only_bacteria = bacteria_cols - archaea_cols
        only_archaea = archaea_cols - bacteria_cols
        
        if only_bacteria:
            print(f"  Columns only in bacteria: {only_bacteria}")
        if only_archaea:
            print(f"  Columns only in archaea: {only_archaea}")
        
        # Find common columns
        common_cols = bacteria_cols & archaea_cols
        print(f"\n  Using {len(common_cols)} common columns")
        
        # Keep only common columns
        bacteria_df = bacteria_df[sorted(common_cols)]
        archaea_df = archaea_df[sorted(common_cols)]
    
    print(f"\nMerging dataframes...")
    # Concatenate bacteria and archaea data
    merged_df = pd.concat([bacteria_df, archaea_df], ignore_index=True)
    
    print(f"  Total assemblies: {len(merged_df)}")
    print(f"  Total columns: {len(merged_df.columns)}")
    
    # Save merged metadata
    print(f"\nSaving merged metadata to {output_file}...")
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    # Show summary statistics
    print("\nMerged metadata statistics:")
    print(f"  Total assemblies: {len(merged_df)}")
    print(f"  Bacteria: {len(bacteria_df)} ({len(bacteria_df)/len(merged_df)*100:.1f}%)")
    print(f"  Archaea: {len(archaea_df)} ({len(archaea_df)/len(merged_df)*100:.1f}%)")
    
    # Show domain distribution from gtdb_taxonomy
    if 'gtdb_taxonomy' in merged_df.columns:
        domains = merged_df['gtdb_taxonomy'].str.split(';').str[0].value_counts()
        print("\nDomain distribution:")
        for domain, count in domains.items():
            print(f"  {domain}: {count}")
    
    print(f"\n✓ Merged metadata saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Merge bacteria and archaea GTDB filtered metadata files'
    )
    parser.add_argument(
        '--bacteria',
        required=True,
        help='Path to bacteria filtered metadata TSV file'
    )
    parser.add_argument(
        '--archaea',
        required=True,
        help='Path to archaea filtered metadata TSV file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to save merged metadata TSV file'
    )
    
    args = parser.parse_args()
    
    merge_metadata(args.bacteria, args.archaea, args.output)


if __name__ == '__main__':
    main()
