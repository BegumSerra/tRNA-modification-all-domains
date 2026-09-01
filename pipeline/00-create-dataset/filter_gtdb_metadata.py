#!/usr/bin/env python3
"""
Filter GTDB metadata to only include assemblies present in our dataset.
This creates a subset of GTDB metadata matching the downloaded proteins.
"""

import pandas as pd
import argparse
from pathlib import Path
import os


def filter_gtdb_metadata(gtdb_metadata_file, accession_list_file, output_file):
    """
    Filter GTDB metadata to only include assemblies from the accession list.
    
    Parameters:
    -----------
    gtdb_metadata_file : str
        Path to full GTDB metadata file
    accession_list_file : str
        Path to file with list of assembly accessions (one per line)
    output_file : str
        Path to save filtered metadata
    """
    print(f"Loading accession list from {accession_list_file}...")
    with open(accession_list_file, 'r') as f:
        accessions = [line.strip() for line in f if line.strip() and not line.startswith('assembly_accession')]
    
    # Remove version numbers and create both GCF and GCA variants for matching
    base_accessions = set()
    for acc in accessions:
        base_acc = acc.rsplit('.', 1)[0]
        base_accessions.add(base_acc)
        
        # Add the alternate version (GCF <-> GCA)
        if base_acc.startswith('GCF_'):
            base_accessions.add(base_acc.replace('GCF_', 'GCA_'))
        elif base_acc.startswith('GCA_'):
            base_accessions.add(base_acc.replace('GCA_', 'GCF_'))
    
    print(f"  Found {len(accessions)} assemblies in list (expanded to {len(base_accessions)} with GCF/GCA variants)")
    
    print(f"Loading GTDB metadata from {gtdb_metadata_file}...")
    # Read in chunks to handle large file
    chunk_size = 10000
    filtered_chunks = []
    total_rows = 0
    matched_rows = 0
    
    for chunk in pd.read_csv(gtdb_metadata_file, sep='\t', chunksize=chunk_size):
        total_rows += len(chunk)
        
        # Match by assembly accession
        # Check both accession column (RefSeq with RS_/GB_ prefix) and ncbi_genbank_assembly_accession (GenBank)
        mask = pd.Series([False] * len(chunk))
        
        if 'accession' in chunk.columns:
            # Remove RS_/GB_ prefix and version number from accession column
            chunk['base_accession'] = chunk['accession'].str.replace('RS_', '').str.replace('GB_', '').str.rsplit('.', n=1).str[0]
            mask = mask | chunk['base_accession'].isin(base_accessions)
        
        if 'ncbi_genbank_assembly_accession' in chunk.columns:
            # Remove version number from GenBank accession
            chunk['base_genbank'] = chunk['ncbi_genbank_assembly_accession'].str.rsplit('.', n=1).str[0]
            mask = mask | chunk['base_genbank'].isin(base_accessions)
        
        filtered_chunk = chunk[mask]
        matched_rows += len(filtered_chunk)
        
        if len(filtered_chunk) > 0:
            filtered_chunks.append(filtered_chunk)
        
        if total_rows % 100000 == 0:
            print(f"  Processed {total_rows} rows, matched {matched_rows} so far...")
    
    print(f"  Total rows processed: {total_rows}")
    print(f"  Matched rows: {matched_rows}")
    
    if not filtered_chunks:
        print("WARNING: No matching assemblies found!")
        return
    
    # Combine all filtered chunks
    filtered_df = pd.concat(filtered_chunks, ignore_index=True)
    
    # Save filtered metadata
    print(f"Saving filtered metadata to {output_file}...")
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nFiltered metadata statistics:")
    print(f"  Input assemblies: {total_rows}")
    print(f"  Output assemblies: {len(filtered_df)}")
    print(f"  Reduction: {(1 - len(filtered_df)/total_rows)*100:.1f}%")
    
    # Show taxonomy coverage
    if 'gtdb_taxonomy' in filtered_df.columns:
        filled_taxonomy = filtered_df['gtdb_taxonomy'].notna().sum()
        print(f"  With GTDB taxonomy: {filled_taxonomy} ({filled_taxonomy/len(filtered_df)*100:.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description='Filter GTDB metadata to only assemblies in the accession list'
    )
    parser.add_argument(
        '--gtdb-metadata',
        required=True,
        help='Path to full GTDB metadata file'
    )
    parser.add_argument(
        '--accession-list',
        required=True,
        help='Path to file with assembly accession list (one per line)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to save filtered metadata'
    )
    
    args = parser.parse_args()
    
    filter_gtdb_metadata(args.gtdb_metadata, args.accession_list, args.output)
    
    print("\nDone!")


if __name__ == '__main__':
    main()
