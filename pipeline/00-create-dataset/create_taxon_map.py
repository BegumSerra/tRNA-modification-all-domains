#!/usr/bin/env python3
"""
Create taxon_map.tsv file mapping protein IDs to taxonomy IDs
for bacteria and archaea assemblies
"""
import pandas as pd
from pathlib import Path
import argparse
import re


def parse_protein_ids_from_faa(faa_file):
    """Extract protein IDs from a FASTA file"""
    protein_ids = []
    with open(faa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract protein ID (first field after >)
                protein_id = line[1:].split()[0]
                protein_ids.append(protein_id)
    return protein_ids


def get_assembly_taxid_map(metadata_file):
    """Create a mapping of assembly accession to taxid from metadata"""
    df = pd.read_csv(metadata_file, sep='\t')
    
    # Create mapping from both RefSeq accession and GenBank accession
    accession_map = {}
    
    # Map RefSeq accession (accession column)
    for _, row in df.iterrows():
        accession = row['accession']
        taxid = row['ncbi_taxid']
        
        # Add with full version
        accession_map[accession] = taxid
        
        # Also add GCF/GCA variant with same version
        if accession.startswith('GCF_'):
            accession_map[accession.replace('GCF_', 'GCA_')] = taxid
        elif accession.startswith('GCA_'):
            accession_map[accession.replace('GCA_', 'GCF_')] = taxid
        
        # Also add without version number (to handle version mismatches)
        base_no_version = accession.rsplit('.', 1)[0]  # Remove .X version
        accession_map[base_no_version] = taxid
        # And the GCF/GCA variant without version
        if base_no_version.startswith('GCF_'):
            accession_map[base_no_version.replace('GCF_', 'GCA_')] = taxid
        elif base_no_version.startswith('GCA_'):
            accession_map[base_no_version.replace('GCA_', 'GCF_')] = taxid
        
        # Also map GenBank accession if available
        if 'ncbi_genbank_assembly_accession' in df.columns and pd.notna(row['ncbi_genbank_assembly_accession']):
            gb_acc = row['ncbi_genbank_assembly_accession']
            accession_map[gb_acc] = taxid
            # Add GCF/GCA variant
            if gb_acc.startswith('GCF_'):
                accession_map[gb_acc.replace('GCF_', 'GCA_')] = taxid
            elif gb_acc.startswith('GCA_'):
                accession_map[gb_acc.replace('GCA_', 'GCF_')] = taxid
            # And without version
            gb_base = gb_acc.rsplit('.', 1)[0]
            accession_map[gb_base] = taxid
            if gb_base.startswith('GCF_'):
                accession_map[gb_base.replace('GCF_', 'GCA_')] = taxid
            elif gb_base.startswith('GCA_'):
                accession_map[gb_base.replace('GCA_', 'GCF_')] = taxid
    
    return accession_map


def create_taxon_map(download_dir, metadata_file, output_file):
    """Create taxon_map.tsv from downloaded protein files and metadata"""
    download_dir = Path(download_dir)
    
    # Get accession -> taxid mapping from metadata
    print(f"Reading metadata from {metadata_file}...")
    accession_to_taxid = get_assembly_taxid_map(metadata_file)
    print(f"Found {len(accession_to_taxid)} assemblies in metadata")
    
    # Find all protein.faa files - try both batch structure and direct structure
    protein_files = list(download_dir.glob("batch_*/ncbi_dataset/data/GC*/protein.faa"))
    if not protein_files:
        # Try direct structure (e.g., ncbi_dataset_130824/data/GC*/protein.faa)
        protein_files = list(download_dir.glob("data/GC*/protein.faa"))
    print(f"\nFound {len(protein_files)} protein files")
    
    if not protein_files:
        print("ERROR: No protein files found!")
        return
    
    # Collect all protein_id -> taxid mappings
    taxon_map = []
    processed_assemblies = 0
    total_proteins = 0
    
    for protein_file in sorted(protein_files):
        # Extract assembly accession from path (e.g., GCF_000023825.1)
        assembly_accession = protein_file.parent.name
        
        # Try exact match first, then without version
        if assembly_accession in accession_to_taxid:
            taxid = accession_to_taxid[assembly_accession]
        else:
            # Try without version
            base_acc = assembly_accession.rsplit('.', 1)[0]
            if base_acc in accession_to_taxid:
                taxid = accession_to_taxid[base_acc]
            else:
                print(f"Warning: {assembly_accession} not found in metadata, skipping...")
                continue
        
        # Parse protein IDs from this file
        protein_ids = parse_protein_ids_from_faa(protein_file)
        
        # Add all proteins from this assembly with the same taxid
        for protein_id in protein_ids:
            taxon_map.append((protein_id, taxid))
        
        processed_assemblies += 1
        total_proteins += len(protein_ids)
        
        if processed_assemblies % 10 == 0:
            print(f"  Processed {processed_assemblies} assemblies, {total_proteins:,} proteins...")
    
    print(f"\nTotal: {processed_assemblies} assemblies, {total_proteins:,} proteins")
    
    # Create DataFrame and save
    df = pd.DataFrame(taxon_map, columns=['protein_id', 'taxid'])
    
    # Save without header
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False, header=False)
    
    print(f"\n✓ Created {output_file}")
    print(f"  Total entries: {len(df):,}")
    print(f"  Unique proteins: {df['protein_id'].nunique():,}")
    print(f"  Unique taxids: {df['taxid'].nunique()}")


def main():
    parser = argparse.ArgumentParser(description='Create taxon_map.tsv from NCBI downloads')
    parser.add_argument('download_dir', help='Directory with downloaded data (contains batch_* folders)')
    parser.add_argument('metadata_file', help='Metadata TSV file with accession and ncbi_taxid columns')
    parser.add_argument('output_file', help='Output taxon_map.tsv file')
    
    args = parser.parse_args()
    
    create_taxon_map(args.download_dir, args.metadata_file, args.output_file)


if __name__ == '__main__':
    main()
