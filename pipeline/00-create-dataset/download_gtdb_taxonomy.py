#!/usr/bin/env python3
"""
Download complete GTDB taxonomy data and merge with existing metadata.

GTDB provides complete taxonomic information for bacterial and archaeal genomes.
This script downloads the latest GTDB metadata and merges it with NCBI metadata.
"""

import pandas as pd
import argparse
import requests
import gzip
import shutil
from pathlib import Path


def download_gtdb_metadata(output_file, domain='archaea'):
    """
    Download GTDB metadata file.
    
    Parameters:
    -----------
    output_file : str
        Path to save the downloaded metadata
    domain : str
        'archaea' or 'bacteria'
    """
    # GTDB Release 226 (latest as of April 2025)
    base_url = "https://data.gtdb.ecogenomic.org/releases/release226/226.0/"
    
    if domain == 'archaea':
        url = f"{base_url}ar53_metadata_r226.tsv.gz"
    elif domain == 'bacteria':
        url = f"{base_url}bac120_metadata_r226.tsv.gz"
    else:
        raise ValueError("Domain must be 'archaea' or 'bacteria'")
    
    print(f"Downloading GTDB {domain} metadata from {url}...")
    
    # Download the file
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    # Save compressed file temporarily
    temp_gz = f"{output_file}.gz"
    with open(temp_gz, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    
    print(f"Decompressing...")
    # Decompress
    with gzip.open(temp_gz, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    # Remove compressed file
    Path(temp_gz).unlink()
    
    print(f"Downloaded GTDB metadata to {output_file}")


def merge_gtdb_with_ncbi(ncbi_metadata, gtdb_metadata, output_file):
    """
    Merge GTDB taxonomy with existing NCBI metadata.
    
    Parameters:
    -----------
    ncbi_metadata : str
        Path to NCBI metadata file
    gtdb_metadata : str
        Path to GTDB metadata file
    output_file : str
        Path to save merged metadata
    """
    print(f"Loading NCBI metadata from {ncbi_metadata}...")
    ncbi_df = pd.read_csv(ncbi_metadata, sep='\t')
    
    print(f"Loading GTDB metadata from {gtdb_metadata}...")
    # Load only necessary columns from GTDB
    gtdb_df = pd.read_csv(gtdb_metadata, sep='\t', usecols=[
        'accession', 'gtdb_taxonomy', 'ncbi_taxonomy', 'ncbi_taxid'
    ])
    
    print(f"Merging metadata...")
    # Merge on accession
    # GTDB uses accessions with version (e.g., GCA_000006155.2)
    # NCBI metadata might have accessions without version
    
    # Try exact match first
    merged = pd.merge(
        ncbi_df,
        gtdb_df[['accession', 'gtdb_taxonomy']],
        on='accession',
        how='left',
        suffixes=('', '_gtdb')
    )
    
    # For rows without match, try matching without version
    no_match = merged[merged['gtdb_taxonomy'].isna()]
    if len(no_match) > 0:
        print(f"  {len(no_match)} accessions without exact match, trying without version...")
        
        # Create accession without version for matching
        gtdb_df['accession_base'] = gtdb_df['accession'].str.rsplit('.', n=1).str[0]
        ncbi_df['accession_base'] = ncbi_df['accession'].str.rsplit('.', n=1).str[0]
        
        # Merge again on base accession
        gtdb_base = gtdb_df[['accession_base', 'gtdb_taxonomy']].drop_duplicates(subset=['accession_base'])
        merged_base = pd.merge(
            ncbi_df[['accession', 'accession_base']],
            gtdb_base,
            on='accession_base',
            how='left'
        )
        
        # Update the rows that were missing
        for idx in no_match.index:
            acc = merged.loc[idx, 'accession']
            gtdb_tax = merged_base[merged_base['accession'] == acc]['gtdb_taxonomy'].values
            if len(gtdb_tax) > 0 and pd.notna(gtdb_tax[0]):
                merged.loc[idx, 'gtdb_taxonomy'] = gtdb_tax[0]
    
    # Update the gtdb_taxonomy column in original dataframe
    if 'gtdb_taxonomy' in merged.columns:
        ncbi_df['gtdb_taxonomy'] = merged['gtdb_taxonomy']
    else:
        ncbi_df = merged
    
    print(f"Saving merged metadata to {output_file}...")
    ncbi_df.to_csv(output_file, sep='\t', index=False)
    
    # Print statistics
    total = len(ncbi_df)
    with_gtdb = ncbi_df['gtdb_taxonomy'].notna().sum()
    print(f"\nStatistics:")
    print(f"  Total genomes: {total}")
    print(f"  With GTDB taxonomy: {with_gtdb} ({with_gtdb/total*100:.1f}%)")
    print(f"  Without GTDB taxonomy: {total-with_gtdb} ({(total-with_gtdb)/total*100:.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description='Download GTDB taxonomy metadata'
    )
    parser.add_argument(
        '--domain',
        choices=['archaea', 'bacteria'],
        required=True,
        help='Domain to process'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to output GTDB metadata file'
    )
    
    args = parser.parse_args()
    
    # Download GTDB metadata
    download_gtdb_metadata(args.output, args.domain)
    
    print("\nDone!")


if __name__ == '__main__':
    main()
