#!/usr/bin/env python3
"""
Create a mapping file with protein IDs, genome IDs, taxon IDs, and species names.
Extracts information from protein.faa files in bacteria and archaea directories.
"""

import os
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd
import sys

def parse_fasta_header(header):
    """
    Parse FASTA header to extract protein ID.
    
    Example headers:
    >WP_052294782.1 MULTISPECIES: tRNA uridine-5-carboxymethylaminomethyl(34) synthesis GTPase MnmE [Acidobacterium]
    >NP_394695.1 queuosine biosynthesis protein QueC [Thermoplasma acidophilum DSM 1728]
    """
    # Extract protein ID (first field after >)
    protein_id_match = re.match(r'>(\S+)', header)
    protein_id = protein_id_match.group(1) if protein_id_match else None
    
    return protein_id

def parse_ncbi_taxonomy(taxonomy_string):
    """
    Parse NCBI taxonomy string to extract individual levels.
    Format: d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
    """
    if pd.isna(taxonomy_string) or not taxonomy_string:
        return {
            'domain': '', 'phylum': '', 'class': '', 
            'order': '', 'family': '', 'genus': '', 'species': ''
        }
    
    taxonomy_dict = {
        'domain': '', 'phylum': '', 'class': '', 
        'order': '', 'family': '', 'genus': '', 'species': ''
    }
    
    parts = taxonomy_string.split(';')
    for part in parts:
        part = part.strip()
        if part.startswith('d__'):
            taxonomy_dict['domain'] = part[3:]
        elif part.startswith('p__'):
            taxonomy_dict['phylum'] = part[3:]
        elif part.startswith('c__'):
            taxonomy_dict['class'] = part[3:]
        elif part.startswith('o__'):
            taxonomy_dict['order'] = part[3:]
        elif part.startswith('f__'):
            taxonomy_dict['family'] = part[3:]
        elif part.startswith('g__'):
            taxonomy_dict['genus'] = part[3:]
        elif part.startswith('s__'):
            taxonomy_dict['species'] = part[3:]
    
    return taxonomy_dict

def load_metadata(metadata_file):
    """
    Load genome metadata from merged_meta_data_bact_arc.tsv.
    Returns a dictionary mapping genome accessions to their metadata.
    Creates multiple lookup keys for GCA/GCF variants and version-less accessions.
    """
    metadata = {}
    
    try:
        df = pd.read_csv(metadata_file, sep='\t')
        
        for _, row in df.iterrows():
            accession = row.get('accession', '')
            taxon_id = row.get('ncbi_taxid', '')
            species_name = row.get('ncbi_organism_name', '')
            ncbi_taxonomy = row.get('ncbi_taxonomy', '')
            
            if pd.notna(accession):
                # Parse NCBI taxonomy into individual levels
                taxonomy_levels = parse_ncbi_taxonomy(ncbi_taxonomy)
                
                meta_info = {
                    'taxon_id': taxon_id if pd.notna(taxon_id) else '',
                    'species_name': species_name if pd.notna(species_name) else '',
                    'domain': taxonomy_levels['domain'],
                    'phylum': taxonomy_levels['phylum'],
                    'class': taxonomy_levels['class'],
                    'order': taxonomy_levels['order'],
                    'family': taxonomy_levels['family'],
                    'genus': taxonomy_levels['genus'],
                    'species': taxonomy_levels['species'],
                    'ncbi_taxonomy': ncbi_taxonomy if pd.notna(ncbi_taxonomy) else ''
                }
                
                # Store with exact accession
                metadata[accession] = meta_info
                
                # Store with base accession (no version)
                base_acc = accession.rsplit('.', 1)[0]
                if base_acc not in metadata:
                    metadata[base_acc] = meta_info
                
                # Store with swapped GCA/GCF
                if accession.startswith('GCA_'):
                    alt_acc = accession.replace('GCA_', 'GCF_')
                    metadata[alt_acc] = meta_info
                    # Also base version of swapped
                    alt_base = alt_acc.rsplit('.', 1)[0]
                    if alt_base not in metadata:
                        metadata[alt_base] = meta_info
                elif accession.startswith('GCF_'):
                    alt_acc = accession.replace('GCF_', 'GCA_')
                    metadata[alt_acc] = meta_info
                    # Also base version of swapped
                    alt_base = alt_acc.rsplit('.', 1)[0]
                    if alt_base not in metadata:
                        metadata[alt_base] = meta_info
        
        print(f"Loaded metadata for {len(df)} genomes (with {len(metadata)} lookup keys)")
    
    except Exception as e:
        print(f"Error loading metadata: {e}")
    
    return metadata

def extract_proteins_from_genome(faa_file, genome_id, metadata):
    """
    Extract all protein IDs from a genome's protein.faa file.
    Returns list of dictionaries with protein info.
    """
    proteins = []
    
    # Get metadata for this genome - try exact match first, then base version, then GCA/GCF swap
    genome_metadata = metadata.get(genome_id, {})
    if not genome_metadata:
        # Try base version (without .X)
        base_id = genome_id.rsplit('.', 1)[0]
        genome_metadata = metadata.get(base_id, {})
    if not genome_metadata:
        # Try swapping GCA/GCF
        if genome_id.startswith('GCA_'):
            alt_id = genome_id.replace('GCA_', 'GCF_')
            genome_metadata = metadata.get(alt_id, {})
            if not genome_metadata:
                # Try base version of swapped
                alt_base = alt_id.rsplit('.', 1)[0]
                genome_metadata = metadata.get(alt_base, {})
        elif genome_id.startswith('GCF_'):
            alt_id = genome_id.replace('GCF_', 'GCA_')
            genome_metadata = metadata.get(alt_id, {})
            if not genome_metadata:
                # Try base version of swapped
                alt_base = alt_id.rsplit('.', 1)[0]
                genome_metadata = metadata.get(alt_base, {})
    
    taxon_id = genome_metadata.get('taxon_id', '')
    species_name = genome_metadata.get('species_name', '')
    domain = genome_metadata.get('domain', '')
    phylum = genome_metadata.get('phylum', '')
    class_name = genome_metadata.get('class', '')
    order = genome_metadata.get('order', '')
    family = genome_metadata.get('family', '')
    genus = genome_metadata.get('genus', '')
    species = genome_metadata.get('species', '')
    ncbi_taxonomy = genome_metadata.get('ncbi_taxonomy', '')
    
    # If species from taxonomy is empty, use species_name (ncbi_organism_name)
    if not species and species_name:
        species = species_name
    
    try:
        with open(faa_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    protein_id = parse_fasta_header(line.strip())
                    
                    if protein_id:
                        proteins.append({
                            'protein_id': protein_id,
                            'genome_id': genome_id,
                            'taxon_id': taxon_id,
                            'domain': domain,
                            'phylum': phylum,
                            'class': class_name,
                            'order': order,
                            'family': family,
                            'genus': genus,
                            'species': species,
                            'ncbi_taxonomy': ncbi_taxonomy
                        })
    
    except Exception as e:
        print(f"Error reading {faa_file}: {e}")
    
    return proteins

def scan_genome_directories(base_dir, metadata):
    """
    Scan through bacteria and archaea directories to extract protein information.
    """
    base_path = Path(base_dir)
    all_proteins = []
    
    # Process bacteria and archaea with their specific dataset paths
    domain_paths = {
        'archaea': base_path / 'archaea' / 'ncbi_dataset_190824' / 'data',
        'bacteria': base_path / 'bacteria' / 'ncbi_dataset_130824' / 'data'
    }
    
    for domain, domain_path in domain_paths.items():
        if not domain_path.exists():
            print(f"Warning: {domain} directory not found at {domain_path}")
            continue
        
        print(f"\nScanning {domain}...")
        
        # Find all protein.faa files
        genome_count = 0
        for faa_file in domain_path.rglob('protein.faa'):
            genome_id = faa_file.parent.name  # e.g., GCF_000007985.2
            
            proteins = extract_proteins_from_genome(faa_file, genome_id, metadata)
            all_proteins.extend(proteins)
            
            genome_count += 1
            if genome_count % 10 == 0:
                print(f"  Processed {genome_count} genomes, {len(all_proteins)} proteins so far...")
        
        print(f"  Completed {domain}: {genome_count} genomes")
    
    return all_proteins

def main():
    # Paths
    base_dir = sys.argv[1]  # e.g., "../../input/database/"
    metadata_file = sys.argv[2]  # e.g., "../../output/table/taxa_dataset/merged_meta_data_bact_arc.tsv"
    output_file = sys.argv[3]  # e.g., "../../output/table/protein_genome_mapping.tsv"
    
    print("Creating protein-genome mapping file...")
    print(f"Base directory: {base_dir}")
    print(f"Metadata file: {metadata_file}")
    print(f"Output file: {output_file}")
    
    # Load metadata
    print("\nLoading genome metadata...")
    metadata = load_metadata(metadata_file)
    
    # Scan genomes and extract protein information
    print("\nScanning genome directories...")
    all_proteins = scan_genome_directories(base_dir, metadata)
    
    if not all_proteins:
        print("\nNo proteins found!")
        return
    
    # Create DataFrame and save
    print(f"\nCreating mapping file with {len(all_proteins)} proteins...")
    df = pd.DataFrame(all_proteins)
    # Create DataFrame and save
    print(f"\nCreating mapping file with {len(all_proteins)} proteins...")
    df = pd.DataFrame(all_proteins)
    
    # Reorder columns
    df = df[['protein_id', 'genome_id', 'taxon_id', 
             'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species',
             'ncbi_taxonomy']]
    
    # Save to TSV
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nDone! Mapping file saved to: {output_file}")
    print(f"\nSummary:")
    print(f"  Total proteins: {len(df):,}")
    print(f"  Unique genomes: {df['genome_id'].nunique():,}")
    print(f"  Unique species: {df['species'].nunique():,}")
    print(f"  Proteins with taxon_id: {df['taxon_id'].notna().sum():,}")
    print(df.head(10).to_string(index=False))

if __name__ == '__main__':
    main()
