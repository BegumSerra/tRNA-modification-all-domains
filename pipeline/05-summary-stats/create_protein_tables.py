#!/usr/bin/env python3
"""
Create comprehensive protein-specific tables combining taxonomy, BLAST results, 
InterProScan domains, and KEGG annotations.

Author: GitHub Copilot
Date: 2025-11-27
"""

import pandas as pd
import argparse
import re
import os
from pathlib import Path
from collections import defaultdict


def load_taxonomy_archaea(metadata_file):
    """Load archaea taxonomy from NCBI taxonomy in GTDB metadata."""
    df = pd.read_csv(metadata_file, sep='\t', usecols=[
        'accession', 'ncbi_taxid', 'ncbi_taxonomy', 'ncbi_organism_name'
    ])
    
    # Note: Accessions are already cleaned by filter_gtdb_metadata.py
    # No need to call clean_accession here
    
    # Parse GTDB taxonomy string (format: d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species)
    def parse_gtdb_taxonomy(tax_string):
        if pd.isna(tax_string) or tax_string == 'none':
            return {'domain': '', 'phylum': '', 'class': '', 'order': '', 'family': '', 'genus': '', 'species': ''}
        
        parts = tax_string.split(';')
        taxonomy = {}
        for part in parts:
            if '__' in part:
                prefix, value = part.split('__', 1)
                if prefix == 'd':
                    taxonomy['domain'] = value
                elif prefix == 'p':
                    taxonomy['phylum'] = value
                elif prefix == 'c':
                    taxonomy['class'] = value
                elif prefix == 'o':
                    taxonomy['order'] = value
                elif prefix == 'f':
                    taxonomy['family'] = value
                elif prefix == 'g':
                    taxonomy['genus'] = value
                elif prefix == 's':
                    taxonomy['species'] = value
        return taxonomy
    
    # Apply parsing to each row
    taxonomy_data = df['ncbi_taxonomy'].apply(parse_gtdb_taxonomy)
    taxonomy_df = pd.DataFrame(taxonomy_data.tolist())
    
    # Combine with accession and taxid
    result = pd.concat([df[['accession', 'ncbi_taxid']], taxonomy_df], axis=1)
    
    # Use NCBI organism name as fallback when GTDB species is empty
    result['species'] = result['species'].replace('', pd.NA)
    result['species'] = result['species'].fillna(df['ncbi_organism_name'])
    
    return result


def load_taxonomy_bacteria(metadata_file):
    """Load bacteria taxonomy from NCBI taxonomy in GTDB metadata."""
    df = pd.read_csv(metadata_file, sep='\t', usecols=[
        'accession', 'ncbi_taxid', 'ncbi_taxonomy', 'ncbi_organism_name'
    ])
    
    # Note: Accessions are already cleaned by filter_gtdb_metadata.py
    # No need to call clean_accession here
    
    # Parse GTDB taxonomy string (format: d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species)
    def parse_gtdb_taxonomy(tax_string):
        if pd.isna(tax_string) or tax_string == 'none':
            return {'domain': '', 'phylum': '', 'class': '', 'order': '', 'family': '', 'genus': '', 'species': ''}
        
        parts = tax_string.split(';')
        taxonomy = {}
        for part in parts:
            if '__' in part:
                prefix, value = part.split('__', 1)
                if prefix == 'd':
                    taxonomy['domain'] = value
                elif prefix == 'p':
                    taxonomy['phylum'] = value
                elif prefix == 'c':
                    taxonomy['class'] = value
                elif prefix == 'o':
                    taxonomy['order'] = value
                elif prefix == 'f':
                    taxonomy['family'] = value
                elif prefix == 'g':
                    taxonomy['genus'] = value
                elif prefix == 's':
                    taxonomy['species'] = value
        return taxonomy
    
    # Apply parsing to each row
    taxonomy_data = df['ncbi_taxonomy'].apply(parse_gtdb_taxonomy)
    taxonomy_df = pd.DataFrame(taxonomy_data.tolist())
    
    # Combine with accession and taxid
    result = pd.concat([df[['accession', 'ncbi_taxid']], taxonomy_df], axis=1)
    
    # Use NCBI organism name as fallback when GTDB species is empty
    result['species'] = result['species'].replace('', pd.NA)
    result['species'] = result['species'].fillna(df['ncbi_organism_name'])
    
    return result


def load_taxonomy_eukarya(taxa_file):
    """Load eukarya taxonomy."""
    df = pd.read_csv(taxa_file, sep='\t')
    
    # Handle different column names - 'acc' or 'id' can be used as protein_id
    if 'acc' in df.columns:
        df = df.rename(columns={'acc': 'protein_id'})
    elif 'id' in df.columns:
        df = df.rename(columns={'id': 'protein_id'})
    
    # Set domain to Eukarya if not present or overwrite
    df['domain'] = 'Eukarya'
    
    # Add missing columns with 'False' if not present
    for col in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
        if col not in df.columns:
            df[col] = 'False'
    
    return df[['protein_id', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']]


def load_taxon_map(taxon_map_file):
    """Load taxon map file (protein_id -> taxid)."""
    if not os.path.exists(taxon_map_file):
        return pd.DataFrame()
    df = pd.read_csv(taxon_map_file, sep='\t', header=None, names=['protein_id', 'taxid'])
    return df


def load_accession_order(accession_file):
    """Load accession list to preserve original order."""
    with open(accession_file, 'r') as f:
        accessions = [line.strip() for line in f if line.strip() and not line.startswith('assembly_accession')]
    # Remove version numbers and create list preserving order
    # Also track both GCF and GCA variants with same order index
    order_list = []
    for acc in accessions:
        base_acc = acc.rsplit('.', 1)[0]
        order_list.append(base_acc)
    return order_list


def load_all_taxonomy(arc_metadata, bact_metadata, euk_taxa, arc_taxon_map, bact_taxon_map, arc_accession_file=None, bact_accession_file=None):
    """Load and combine all taxonomy data with taxon maps."""
    print("Loading taxonomy data...")
    arc_tax = load_taxonomy_archaea(arc_metadata)
    print(f"  Loaded {len(arc_tax)} archaea taxonomy entries")
    
    bact_tax = load_taxonomy_bacteria(bact_metadata)
    print(f"  Loaded {len(bact_tax)} bacteria taxonomy entries")
    
    euk_tax = load_taxonomy_eukarya(euk_taxa)
    print(f"  Loaded {len(euk_tax)} eukarya taxonomy entries")
    
    # Load taxon maps (protein_id -> taxid)
    arc_map = load_taxon_map(arc_taxon_map)
    print(f"  Loaded {len(arc_map)} archaea protein->taxid mappings")
    
    bact_map = load_taxon_map(bact_taxon_map)
    print(f"  Loaded {len(bact_map)} bacteria protein->taxid mappings")
    
    # Merge taxon maps with taxonomy for archaea
    # taxon_map has: protein_id -> taxid
    # metadata has: ncbi_taxid -> taxonomy info
    # Use 'right' join to keep ALL species from metadata, even without proteins
    if not arc_map.empty and 'ncbi_taxid' in arc_tax.columns:
        arc_tax['ncbi_taxid'] = arc_tax['ncbi_taxid'].astype(str)
        arc_map['taxid'] = arc_map['taxid'].astype(str)
        arc_with_proteins = pd.merge(arc_map, arc_tax, left_on='taxid', right_on='ncbi_taxid', how='right')
        # For species without proteins, fill protein_id with a dummy value
        arc_with_proteins['protein_id'] = arc_with_proteins['protein_id'].fillna('NO_PROTEIN_' + arc_with_proteins['ncbi_taxid'])
        print(f"  Merged archaea: {len(arc_with_proteins)} protein-taxonomy mappings")
    else:
        arc_with_proteins = pd.DataFrame()
        print("  Warning: Could not merge archaea data")
    
    # Merge taxon maps with taxonomy for bacteria
    # Use 'right' join to keep ALL species from metadata, even without proteins
    if not bact_map.empty and 'ncbi_taxid' in bact_tax.columns:
        bact_tax['ncbi_taxid'] = bact_tax['ncbi_taxid'].astype(str)
        bact_map['taxid'] = bact_map['taxid'].astype(str)
        bact_with_proteins = pd.merge(bact_map, bact_tax, left_on='taxid', right_on='ncbi_taxid', how='right')
        # For species without proteins, fill protein_id with a dummy value
        bact_with_proteins['protein_id'] = bact_with_proteins['protein_id'].fillna('NO_PROTEIN_' + bact_with_proteins['ncbi_taxid'])
        print(f"  Merged bacteria: {len(bact_with_proteins)} protein-taxonomy mappings")
    else:
        bact_with_proteins = pd.DataFrame()
        print("  Warning: Could not merge bacteria data")
    
    # For eukarya, protein_id column already contains the EP* identifiers
    # No additional mapping needed
    
    # Load accession order if provided
    arc_order = load_accession_order(arc_accession_file) if arc_accession_file else []
    bact_order = load_accession_order(bact_accession_file) if bact_accession_file else []
    
    # Create order mappings - also map GCF<->GCA variants to same index
    # Archaea indices start at 0
    arc_order_map = {}
    for idx, acc in enumerate(arc_order):
        arc_order_map[acc] = idx
        # Add alternate version (GCF <-> GCA)
        if acc.startswith('GCF_'):
            arc_order_map[acc.replace('GCF_', 'GCA_')] = idx
        elif acc.startswith('GCA_'):
            arc_order_map[acc.replace('GCA_', 'GCF_')] = idx
    
    # Bacteria indices start after archaea (offset by length of arc_order)
    bact_order_map = {}
    bact_offset = len(arc_order)
    for idx, acc in enumerate(bact_order):
        bact_order_map[acc] = idx + bact_offset
        # Add alternate version (GCF <-> GCA)
        if acc.startswith('GCF_'):
            bact_order_map[acc.replace('GCF_', 'GCA_')] = idx + bact_offset
        elif acc.startswith('GCA_'):
            bact_order_map[acc.replace('GCA_', 'GCF_')] = idx + bact_offset
    
    # Add accession column for ordering by merging with taxonomy
    if not arc_with_proteins.empty and arc_order and 'accession' in arc_tax.columns:
        # Create temporary mapping from taxid to accession
        taxid_to_acc = dict(zip(arc_tax['ncbi_taxid'].astype(str), arc_tax['accession']))
        # Use ncbi_taxid (not taxid) since taxid is NaN for species without proteins
        arc_with_proteins['accession'] = arc_with_proteins['ncbi_taxid'].map(taxid_to_acc)
        arc_with_proteins['sort_order'] = arc_with_proteins['accession'].apply(lambda x: arc_order_map.get(x.rsplit('.', 1)[0], 999999) if pd.notna(x) else 999999)
    
    if not bact_with_proteins.empty and bact_order and 'accession' in bact_tax.columns:
        # Create temporary mapping from taxid to accession
        taxid_to_acc = dict(zip(bact_tax['ncbi_taxid'].astype(str), bact_tax['accession']))
        # Use ncbi_taxid (not taxid) since taxid is NaN for species without proteins
        bact_with_proteins['accession'] = bact_with_proteins['ncbi_taxid'].map(taxid_to_acc)
        bact_with_proteins['sort_order'] = bact_with_proteins['accession'].apply(lambda x: bact_order_map.get(x.rsplit('.', 1)[0], 999999) if pd.notna(x) else 999999)
    
    # Combine all taxonomy data with protein mappings
    taxonomy_list = []
    if not arc_with_proteins.empty:
        cols = ['protein_id', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        if 'sort_order' in arc_with_proteins.columns:
            cols.append('sort_order')
        taxonomy_list.append(arc_with_proteins[cols])
    if not bact_with_proteins.empty:
        cols = ['protein_id', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        if 'sort_order' in bact_with_proteins.columns:
            cols.append('sort_order')
        taxonomy_list.append(bact_with_proteins[cols])
    # For eukarya, assign sort_order based on taxonomic grouping (after archaea and bacteria)
    euk_with_order = euk_tax[['protein_id', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']].copy()
    
    # Create sort_order for eukarya based on taxonomic hierarchy
    # Start after archaea (90) and bacteria (110) = 200
    euk_offset = len(arc_order) + len(bact_order)
    
    # Group eukarya by taxonomy and assign sequential order
    euk_grouped = euk_with_order.groupby(['phylum', 'class', 'order', 'family', 'genus', 'species']).size().reset_index(name='count')
    # Sort by taxonomic hierarchy
    euk_grouped = euk_grouped.sort_values(['phylum', 'class', 'order', 'family', 'genus', 'species'])
    # Assign order numbers
    euk_order_map = {}
    for idx, row in euk_grouped.iterrows():
        key = (row['phylum'], row['class'], row['order'], row['family'], row['genus'], row['species'])
        euk_order_map[key] = idx + euk_offset
    
    # Map sort_order back to euk_with_order
    euk_with_order['sort_order'] = euk_with_order.apply(
        lambda row: euk_order_map.get((row['phylum'], row['class'], row['order'], row['family'], row['genus'], row['species']), 999999),
        axis=1
    )
    
    taxonomy_list.append(euk_with_order)
    
    all_taxonomy = pd.concat(taxonomy_list, ignore_index=True)
    print(f"  Total combined taxonomy mappings: {len(all_taxonomy)}")
    
    return all_taxonomy
    
    return all_taxonomy


def extract_protein_id(blast_hit):
    """
    Extract protein ID from BLAST hit identifier.
    The protein ID is typically the part between pipes: database|PROTEIN_ID|...
    Examples:
    - ref|WP_123456.1| -> WP_123456.1
    - gb|UCE11671.1| -> UCE11671.1
    - sp|Q7Z6V5|ADAT2_HUMAN -> Q7Z6V5
    - plain_id -> plain_id
    """
    if '|' not in blast_hit:
        return blast_hit
    
    parts = blast_hit.split('|')
    
    # If we have at least 2 parts (e.g., "db|ID" or "db|ID|..."), 
    # the protein ID is typically at index 1 (second part)
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    
    # Fallback: return last non-empty part
    for part in reversed(parts):
        if part:
            return part
    
    return blast_hit


def extract_ep_number(protein_id):
    """
    Extract EP number from protein ID.
    For eukarya proteins with format like 'EP00074_Homo_sapiens_P035512', 
    extract just 'EP00074' to match with euk_taxa.tsv
    """
    if protein_id.startswith('EP'):
        # Split by underscore and take the first part (EP#####)
        return protein_id.split('_')[0]
    return protein_id


def load_blast_results(blast_file, protein_name):
    """Load BLAST results for a specific protein."""
    if not os.path.exists(blast_file):
        return pd.DataFrame()
    
    # BLAST output format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids
    blast_df = pd.read_csv(blast_file, sep='\t', header=None, names=[
        'query', 'subject', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'taxid'
    ])
    
    # Extract protein ID from subject
    blast_df['protein_id'] = blast_df['subject'].apply(extract_protein_id)
    
    # Select relevant columns
    blast_df = blast_df[['protein_id', 'length', 'evalue']].copy()
    blast_df = blast_df.rename(columns={'protein_id': 'id', 'evalue': 'e-value'})
    
    # Remove duplicates, keep best e-value
    blast_df = blast_df.sort_values('e-value').drop_duplicates(subset=['id'], keep='first')
    
    return blast_df


def load_domain_results(domain_file):
    """Load InterProScan domain results."""
    if not os.path.exists(domain_file):
        return pd.DataFrame()
    
    # Domain file format from InterProScan
    domain_df = pd.read_csv(domain_file, sep='\t', header=None)
    
    # Check number of columns and parse accordingly
    if domain_df.shape[1] >= 13:
        domain_df.columns = [
            'id', 'md5', 'length', 'analysis', 'signature_accession', 
            'signature_description', 'start', 'stop', 'score', 'status', 
            'date', 'interpro_accession', 'interpro_description'
        ] + [f'col_{i}' for i in range(13, domain_df.shape[1])]
    else:
        return pd.DataFrame()
    
    # Get all unique protein IDs
    all_protein_ids = domain_df['id'].unique()
    
    # Filter to only Pfam domains (starting with PF)
    pfam_df = domain_df[domain_df['signature_accession'].str.startswith('PF', na=False)].copy()
    
    # Combine domain information
    pfam_df['domain'] = pfam_df.apply(
        lambda row: f"{row['signature_accession']}:{row['signature_description']}" 
        if pd.notna(row['signature_description']) else row['signature_accession'],
        axis=1
    )
    
    # Group by ID and concatenate domains
    domain_grouped = pfam_df.groupby('id')['domain'].apply(lambda x: '; '.join(x.unique())).reset_index()
    
    # Create DataFrame with all protein IDs
    all_ids_df = pd.DataFrame({'id': all_protein_ids})
    
    # Merge to include all proteins, marking those without PFAM domains as "False"
    result = pd.merge(all_ids_df, domain_grouped, on='id', how='left')
    result['domain'] = result['domain'].fillna('False')
    
    return result


def load_kegg_annotations(kegg_file):
    """Load KEGG annotations for a protein."""
    if not os.path.exists(kegg_file):
        return pd.DataFrame()
    
    # KEGG file format: protein_id\tKEGG_id
    kegg_df = pd.read_csv(kegg_file, sep='\t', header=None, names=['id', 'kegg'])
    
    return kegg_df


def create_protein_table(protein_name, taxonomy_df, blast_dirs, domain_dirs, kegg_file, output_file, fasta_file=None):
    """
    Create a comprehensive table for a specific protein with all species listed.
    
    Parameters:
    -----------
    protein_name : str
        Name of the protein (e.g., 'ArcS', 'Elp3')
    taxonomy_df : pd.DataFrame
        Combined taxonomy data with protein_id mappings
    blast_dirs : dict
        Dictionary with keys 'arc', 'bact', 'euk' and values as directory paths
    domain_dirs : dict
        Dictionary with keys 'arc', 'bact', 'euk' and values as directory paths
    kegg_file : str
        Path to KEGG annotation file for this protein
    output_file : str
        Path to output table file
    fasta_file : str, optional
        Path to FASTA file containing sequences to include in table
    """
    
    # Start with all UNIQUE species from taxonomy
    # Get unique species by dropping duplicate taxonomy rows, keeping MINIMUM sort_order for each species
    cols_to_keep = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    if 'sort_order' in taxonomy_df.columns:
        cols_to_keep.append('sort_order')
        # Group by taxonomy columns and take minimum sort_order for each unique species
        base_table = taxonomy_df[cols_to_keep].groupby(
            ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'], 
            as_index=False
        )['sort_order'].min()
    else:
        base_table = taxonomy_df[cols_to_keep].drop_duplicates(subset=['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    
    # Filter out rows with empty/False taxonomy (invalid entries)
    base_table = base_table[base_table['domain'] != '']
    base_table = base_table[base_table['domain'] != 'False']
    base_table = base_table[base_table['species'] != '']
    base_table = base_table[base_table['species'] != 'False']
    
    # Rename domain to domain_tax
    base_table = base_table.rename(columns={'domain': 'domain_tax'})
    
    # Sort by sort_order if available, otherwise by taxonomy
    if 'sort_order' in base_table.columns:
        base_table = base_table.sort_values('sort_order')
        base_table = base_table.drop(columns=['sort_order'])
    else:
        # Sort by domain to get Archaea, then Bacteria, then Eukarya
        domain_order = {'Archaea': 0, 'Bacteria': 1, 'Eukarya': 2}
        base_table['domain_sort'] = base_table['domain_tax'].map(domain_order)
        base_table = base_table.sort_values(['domain_sort', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        base_table = base_table.drop(columns=['domain_sort'])
    
    # Reset index
    base_table = base_table.reset_index(drop=True)
    
    # Initialize columns for protein data with False
    base_table['id'] = 'False'
    base_table['length'] = 'False'
    base_table['e-value'] = 'False'
    base_table['domain'] = 'False'
    base_table['kegg'] = 'False'
    
    all_data = []
    
    # Process BLAST results from cross-domain searches
    # BLAST files are organized as query-target (e.g., euk-arc, arc-bact, etc.)
    # We need to check all possible combinations to get hits in all target domains
    search_patterns = [
        ('arc', 'arc'), ('arc', 'bact'), ('arc', 'euk'),
        ('bact', 'arc'), ('bact', 'bact'), ('bact', 'euk'),
        ('euk', 'arc'), ('euk', 'bact'), ('euk', 'euk')
    ]
    
    # Store merged directory path for later use
    merged_blast_dir = os.path.join(os.path.dirname(blast_dirs['arc']), 'merged')
    merged_blast_file = os.path.join(merged_blast_dir, f"{protein_name}.tab")
    
    # Add merged pattern if the file exists
    if os.path.exists(merged_blast_file):
        search_patterns.append(('merged', 'merged'))
    
    for query_domain, target_domain in search_patterns:
        # Handle merged directory specially
        if query_domain == 'merged':
            blast_file = merged_blast_file
        else:
            # Try the cross-domain directory pattern (e.g., euk-arc/)
            cross_domain_dir = os.path.join(os.path.dirname(blast_dirs['arc']), f"{query_domain}-{target_domain}")
            blast_file = os.path.join(cross_domain_dir, f"{protein_name}_unique_filtered.tab")
            
            if not os.path.exists(blast_file):
                # Try the simple domain directory pattern (e.g., arc/)
                blast_file = os.path.join(blast_dirs[target_domain], f"{protein_name}_unique_filtered.tab")
            
            if not os.path.exists(blast_file):
                continue
            
        blast_df = load_blast_results(blast_file, protein_name)
        
        if blast_df.empty:
            continue
        
        # Load domain results - check multiple possible locations
        domain_df = pd.DataFrame()
        
        # Handle merged case specially
        if query_domain == 'merged':
            merged_domain_dir = os.path.join(os.path.dirname(domain_dirs['arc']), 'merged')
            domain_file = os.path.join(merged_domain_dir, f"{protein_name}.tsv")
            if os.path.exists(domain_file):
                domain_df = load_domain_results(domain_file)
        else:
            # Try cross-domain directory first
            cross_domain_dir_domains = os.path.join(os.path.dirname(domain_dirs['arc']), f"{query_domain}-{target_domain}")
            domain_file = os.path.join(cross_domain_dir_domains, f"{protein_name}_filtered_length.tsv")
            
            if os.path.exists(domain_file):
                domain_df = load_domain_results(domain_file)
            
            # If not found, try the simple domain directory pattern
            if domain_df.empty:
                domain_file = os.path.join(domain_dirs[target_domain], f"{protein_name}_filtered_length.tsv")
                if os.path.exists(domain_file):
                    domain_df = load_domain_results(domain_file)
            
            # If still not found, try merged directory
            if domain_df.empty:
                merged_domain_dir = os.path.join(os.path.dirname(domain_dirs['arc']), 'merged')
                domain_file = os.path.join(merged_domain_dir, f"{protein_name}.tsv")
                if os.path.exists(domain_file):
                    domain_df = load_domain_results(domain_file)
        
        # Merge BLAST and domain data
        if not domain_df.empty:
            data = pd.merge(blast_df, domain_df, on='id', how='left')
        else:
            data = blast_df.copy()
            data['domain'] = ''
        
        all_data.append(data)
    
    # Combine all protein data
    if not all_data:
        print(f"No BLAST data found for protein {protein_name}, creating table with all False values")
        combined_data = pd.DataFrame()
    else:
        combined_data = pd.concat(all_data, ignore_index=True)
    
    # Filter based on FASTA file if provided
    if fasta_file and os.path.exists(fasta_file):
        fasta_ids = set()
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Extract ID from FASTA header (first word after >)
                    fasta_ids.add(line.split()[0][1:])
        
        if fasta_ids and not combined_data.empty:
            # Filter combined_data to only include IDs in FASTA
            before_filter = len(combined_data)
            combined_data = combined_data[combined_data['id'].isin(fasta_ids)]
            # Keep only first occurrence of each ID (remove duplicates)
            combined_data = combined_data.drop_duplicates(subset=['id'], keep='first')
            after_filter = len(combined_data)
            print(f"  Filtered table to match FASTA: {before_filter} rows -> {after_filter} rows ({len(fasta_ids)} IDs in FASTA)")
    
    # If we have protein data, merge it with the base table
    if not combined_data.empty:
        # Load KEGG annotations
        kegg_df = load_kegg_annotations(kegg_file)
        
        # Merge with KEGG data
        if not kegg_df.empty:
            combined_data = pd.merge(combined_data, kegg_df, on='id', how='left')
        else:
            combined_data['kegg'] = ''
        
        # Fill NaN values with empty strings
        combined_data = combined_data.fillna('')
        
        # Remove duplicates with _1 suffix when matching entry without suffix exists
        duplicate_ids = []
        ids_with_suffix = combined_data[combined_data['id'].str.contains(r'_\d+$', regex=True)]['id'].tolist()
        
        for id_with_suffix in ids_with_suffix:
            base_id = re.sub(r'_\d+$', '', id_with_suffix)
            if base_id in combined_data['id'].values:
                suffix_entry = combined_data[combined_data['id'] == id_with_suffix]
                # Check if we can determine domain from the ID
                has_match = base_id in combined_data['id'].values
                if has_match:
                    duplicate_ids.append(id_with_suffix)
        
        # Log removed duplicates
        if duplicate_ids:
            print(f"\nRemoving {len(duplicate_ids)} duplicate entries with suffix:")
            for dup_id in duplicate_ids:
                print(f"  - {dup_id}")
            combined_data = combined_data[~combined_data['id'].isin(duplicate_ids)]
        
        # Prepare protein data for merging with taxonomy
        # For eukarya proteins (EP*), extract just the EP number for matching
        combined_data['id_for_match'] = combined_data['id'].apply(extract_ep_number)
        
        # Merge protein data with full taxonomy to get species information
        taxonomy_for_merge = taxonomy_df.copy()
        taxonomy_for_merge['id_for_match'] = taxonomy_for_merge['protein_id'].apply(extract_ep_number)
        
        protein_with_tax = pd.merge(
            combined_data,
            taxonomy_for_merge[['id_for_match', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']],
            on='id_for_match',
            how='left',
            suffixes=('_pfam', '_tax')
        )
        
        # Rename columns: domain_tax is taxonomic domain, domain_pfam is Pfam domain
        if 'domain_pfam' in protein_with_tax.columns:
            protein_with_tax = protein_with_tax.rename(columns={'domain_pfam': 'domain'})
        if 'domain_tax' in protein_with_tax.columns:
            pass  # Already named correctly
        elif 'domain' in protein_with_tax.columns and 'domain_pfam' not in protein_with_tax.columns:
            # If no collision, 'domain' is from taxonomy, rename it
            protein_with_tax = protein_with_tax.rename(columns={'domain': 'domain_tax'})
        
        # For each species in base_table, find if there's matching protein data
        result_rows = []
        for idx, species_row in base_table.iterrows():
            # Skip rows with False or empty values in taxonomy
            if species_row['domain_tax'] == 'False' or species_row['domain_tax'] == '' or pd.isna(species_row['domain_tax']):
                continue
            
            # Find matching proteins for this species
            matching_proteins = protein_with_tax[
                (protein_with_tax['domain_tax'] == species_row['domain_tax']) &
                (protein_with_tax['phylum'] == species_row['phylum']) &
                (protein_with_tax['class'] == species_row['class']) &
                (protein_with_tax['order'] == species_row['order']) &
                (protein_with_tax['family'] == species_row['family']) &
                (protein_with_tax['genus'] == species_row['genus']) &
                (protein_with_tax['species'] == species_row['species'])
            ]
            
            if not matching_proteins.empty:
                # Keep all matching proteins (changed from taking only first)
                for _, match in matching_proteins.iterrows():
                    result_row = {
                        'domain_tax': species_row['domain_tax'],
                        'phylum': species_row['phylum'],
                        'class': species_row['class'],
                        'order': species_row['order'],
                        'family': species_row['family'],
                        'genus': species_row['genus'],
                        'species': species_row['species'],
                        'id': match['id'],
                        'length': match['length'],
                        'e-value': match['e-value'],
                        'domain': match['domain'],
                        'kegg': match['kegg']
                    }
                    result_rows.append(result_row)
            else:
                # No protein found for this species, use False
                result_row = {
                    'domain_tax': species_row['domain_tax'],
                    'phylum': species_row['phylum'],
                    'class': species_row['class'],
                    'order': species_row['order'],
                    'family': species_row['family'],
                    'genus': species_row['genus'],
                    'species': species_row['species'],
                    'id': 'False',
                    'length': 'False',
                    'e-value': 'False',
                    'domain': 'False',
                    'kegg': 'False'
                }
                result_rows.append(result_row)
        
        result_table = pd.DataFrame(result_rows)
    else:
        # No protein data found, use base table with all False
        result_table = base_table.copy()
    
    # Replace empty strings with 'False'
    result_table = result_table.fillna('False')
    for col in result_table.columns:
        result_table[col] = result_table[col].replace('', 'False')
    
    # Reorder columns - only keep domain, phylum, order, family, species, and protein data
    column_order = ['domain_tax', 'phylum', 'order', 'family', 'species',
                   'id', 'length', 'e-value', 'domain', 'kegg']
    
    result_table = result_table[column_order]
    
    # Remove duplicate protein IDs (can occur when same protein ID maps to multiple species)
    # Keep only first occurrence of each ID while maintaining original order
    if 'id' in result_table.columns and not result_table.empty:
        # Add original index to maintain order
        result_table = result_table.reset_index(drop=True)
        result_table['_orig_idx'] = result_table.index
        
        # For non-False IDs, mark duplicates (keep first)
        seen_ids = set()
        rows_to_keep = []
        for idx, row in result_table.iterrows():
            if row['id'] == 'False':
                # Always keep False rows (no protein found)
                rows_to_keep.append(idx)
            elif row['id'] not in seen_ids:
                # First occurrence of this protein ID
                seen_ids.add(row['id'])
                rows_to_keep.append(idx)
            # else: duplicate, skip
        
        result_table = result_table.loc[rows_to_keep].drop(columns=['_orig_idx'])
    
    # Save to file
    result_table.to_csv(output_file, sep='\t', index=False)
    print(f"Created table for {protein_name}: {output_file}")
    print(f"  Total species: {len(result_table)}")
    print(f"  Species with protein: {(result_table['length'] != 'False').sum()}")
    print(f"  With domains: {(result_table['domain'] != 'False').sum()}")
    print(f"  With KEGG: {(result_table['kegg'] != 'False').sum()}")


def get_protein_list(protein_kegg_file):
    """Get list of all proteins from protein_kegg.tsv."""
    df = pd.read_csv(protein_kegg_file)
    return df['protein'].unique().tolist()


def main():
    parser = argparse.ArgumentParser(
        description='Create protein-specific tables with taxonomy and annotations'
    )
    parser.add_argument('--protein', type=str, help='Specific protein name (e.g., ArcS, Elp3)')
    parser.add_argument('--all', action='store_true', help='Process all proteins')
    parser.add_argument('--arc-metadata', type=str, required=True,
                       help='Path to archaea GTDB metadata file')
    parser.add_argument('--bact-metadata', type=str, required=True,
                       help='Path to bacteria GTDB metadata file')
    parser.add_argument('--euk-taxa', type=str, required=True,
                       help='Path to eukarya taxonomy file')
    parser.add_argument('--arc-taxon-map', type=str, required=True,
                       help='Path to archaea taxon map file')
    parser.add_argument('--bact-taxon-map', type=str, required=True,
                       help='Path to bacteria taxon map file')
    parser.add_argument('--arc-accessions', type=str, required=False,
                       help='Path to archaea accession list file (for ordering)')
    parser.add_argument('--bact-accessions', type=str, required=False,
                       help='Path to bacteria accession list file (for ordering)')
    parser.add_argument('--blast-base', type=str, required=True,
                       help='Base directory for BLAST results')
    parser.add_argument('--domain-base', type=str, required=True,
                       help='Base directory for domain results')
    parser.add_argument('--kegg-dir', type=str, required=True,
                       help='Directory containing KEGG annotation files')
    parser.add_argument('--protein-list', type=str, required=True,
                       help='File containing list of all proteins')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory for protein tables')
    parser.add_argument('--fasta-dir', type=str, required=False,
                       help='Directory containing FASTA files (e.g., output/seqs/homologs)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load taxonomy data
    print("Loading taxonomy data...")
    taxonomy_df = load_all_taxonomy(
        args.arc_metadata, 
        args.bact_metadata, 
        args.euk_taxa,
        args.arc_taxon_map,
        args.bact_taxon_map,
        args.arc_accessions,
        args.bact_accessions
    )
    print(f"Loaded taxonomy mappings for {len(taxonomy_df)} proteins")
    
    # Define BLAST and domain directories for each domain
    blast_dirs = {
        'arc': os.path.join(args.blast_base, 'arc-arc'),
        'bact': os.path.join(args.blast_base, 'bact-bact'),
        'euk': os.path.join(args.blast_base, 'euk-euk')
    }
    
    domain_dirs = {
        'arc': os.path.join(args.domain_base, 'arc-arc'),
        'bact': os.path.join(args.domain_base, 'bact-bact'),
        'euk': os.path.join(args.domain_base, 'euk-euk')
    }
    
    # Determine which proteins to process
    if args.all:
        proteins = get_protein_list(args.protein_list)
        print(f"Processing all {len(proteins)} proteins")
    elif args.protein:
        proteins = [args.protein]
        print(f"Processing protein: {args.protein}")
    else:
        parser.error("Must specify either --protein or --all")
    
    # Process each protein
    for protein in proteins:
        print(f"\nProcessing {protein}...")
        kegg_file = os.path.join(args.kegg_dir, f"{protein}_kegg.txt")
        output_file = os.path.join(args.output_dir, f"{protein}_table.tsv")
        
        # Determine FASTA file path if fasta_dir provided
        fasta_file = None
        if args.fasta_dir:
            fasta_file = os.path.join(args.fasta_dir, f"{protein}_nogaps.fa")
            if not os.path.exists(fasta_file):
                print(f"  Warning: FASTA file not found: {fasta_file}")
                fasta_file = None
        
        try:
            create_protein_table(
                protein_name=protein,
                taxonomy_df=taxonomy_df,
                blast_dirs=blast_dirs,
                domain_dirs=domain_dirs,
                kegg_file=kegg_file,
                output_file=output_file,
                fasta_file=fasta_file
            )
        except Exception as e:
            import traceback
            print(f"Error processing {protein}: {str(e)}")
            traceback.print_exc()
            continue
    
    print("\nDone!")


if __name__ == '__main__':
    main()
