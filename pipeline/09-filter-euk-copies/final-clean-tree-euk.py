#!/usr/bin/env python3
"""
Read a phylogenetic tree in Nexus format and filter eukaryotic sequences.
When multiple eukaryotic sequences from the same species cluster together,
keep only the longest one.

Eukaryotic IDs are identified by their style (typically contain underscores or specific patterns).
"""

import sys
from Bio import Phylo
from Bio import SeqIO
import re
from collections import defaultdict

def is_eukaryotic_id(seq_id):
    """
    Identify if a sequence ID is from eukaryotes based on naming pattern.
    Eukaryotic IDs are in format: EP00295_Genus_species_P002334
    """
    # Eukaryotic IDs start with EP followed by digits
    if seq_id.startswith('EP') and len(seq_id) > 2 and seq_id[2:7].isdigit():
        return True
    return False

def extract_species_name(seq_id):
    """
    Extract species name from sequence ID.
    Format: EP00295_Genus_species_P002334 -> Genus_species
    """
    # Split by underscore
    parts = seq_id.split('_')
    if len(parts) >= 3 and seq_id.startswith('EP') and parts[0][2:].isdigit():
        # Skip EP##### (first part) and P##### (last part), keep genus and species
        # Handle case where there might be more underscores in species name
        # Take everything between the first EP##### and last P##### part
        if parts[-1].startswith('P') and parts[-1][1:].isdigit():
            # Standard format: EP#####_Genus_species_P#####
            species_parts = parts[1:-1]
            return '_'.join(species_parts)
        else:
            # No P##### at end, take everything after EP#####
            return '_'.join(parts[1:])
    print(f"    WARNING: Could not extract species from ID: {seq_id}")
    return None

def get_sequence_lengths(fasta_file):
    """
    Read fasta file and return dictionary of sequence ID to length.
    """
    lengths = {}
    with open(fasta_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            lengths[record.id] = len(record.seq)
    return lengths

def get_all_descendants(tree, clade):
    """
    Get all terminal (leaf) descendants of a clade.
    """
    terminals = []
    for terminal in clade.get_terminals():
        terminals.append(terminal.name)
    return terminals

def find_euk_clusters(tree):
    """
    Find clusters of eukaryotic sequences from the same species that are phylogenetically close.
    Consider sequences clustered if:
    1. They dominate a clade (>50% of eukaryotic sequences), OR
    2. They're in a tight cluster (≤6 total eukaryotic sequences), OR  
    3. They appear as sister taxa (same species on both sides of binary split), OR
    4. Sister clades both contain sequences from the same species
    Returns dict: species_name -> list of (clade, seq_ids)
    """
    euk_clusters = defaultdict(list)
    processed_sequences = set()  # Track sequences already in clusters to avoid duplicates
    
    # Traverse all internal nodes from bottom up (postorder) to find tight clusters first
    for clade in tree.find_clades(order='postorder'):
        if clade.is_terminal():
            continue
        
        # Get all terminal descendants of this clade
        all_terminals = get_all_descendants(tree, clade)
        
        # Filter for eukaryotic sequences
        euk_descendants = [seq_id for seq_id in all_terminals if is_eukaryotic_id(seq_id)]
        
        if len(euk_descendants) < 2:
            continue
        
        # Skip if any of these sequences are already processed in a smaller cluster
        if any(seq_id in processed_sequences for seq_id in euk_descendants):
            continue
            
        # Group by species
        species_groups = defaultdict(list)
        for seq_id in euk_descendants:
            species = extract_species_name(seq_id)
            if species:  # Make sure species extraction succeeded
                species_groups[species].append(seq_id)
        
        # Find species with multiple sequences in this clade
        for species, seq_ids in species_groups.items():
            if len(seq_ids) >= 2:
                # Check different clustering criteria:
                species_fraction = len(seq_ids) / len(euk_descendants)
                is_dominant = species_fraction > 0.5
                is_tight_cluster = len(euk_descendants) <= 6
                
                # Check for sister taxa pattern
                is_sister_taxa = False
                if len(clade.clades) == 2:
                    left_terminals = get_all_descendants(tree, clade.clades[0])
                    right_terminals = get_all_descendants(tree, clade.clades[1])
                    
                    left_species_seqs = [seq_id for seq_id in seq_ids if seq_id in left_terminals]
                    right_species_seqs = [seq_id for seq_id in seq_ids if seq_id in right_terminals]
                    
                    if len(left_species_seqs) > 0 and len(right_species_seqs) > 0:
                        is_sister_taxa = True
                
                # Check for sister clades containing same species
                is_sister_clades_same_species = False
                if len(clade.clades) == 2 and len(seq_ids) >= 4:
                    # If we have 4+ sequences of same species in a binary split,
                    # and both sides contain this species, treat as one cluster
                    left_terminals = get_all_descendants(tree, clade.clades[0])
                    right_terminals = get_all_descendants(tree, clade.clades[1])
                    
                    left_species_seqs = [seq_id for seq_id in seq_ids if seq_id in left_terminals]
                    right_species_seqs = [seq_id for seq_id in seq_ids if seq_id in right_terminals]
                    
                    if len(left_species_seqs) >= 2 and len(right_species_seqs) >= 2:
                        is_sister_clades_same_species = True
                
                # Consider it a cluster if any criteria is met
                if is_dominant or is_tight_cluster or is_sister_taxa or is_sister_clades_same_species:
                    euk_clusters[species].append((clade, seq_ids))
                    processed_sequences.update(seq_ids)
                    
                    criteria = []
                    if is_dominant: criteria.append(f"dominant({species_fraction:.2f})")
                    if is_tight_cluster: criteria.append(f"tight({len(euk_descendants)})")
                    if is_sister_taxa: criteria.append("sister_taxa")
                    if is_sister_clades_same_species: criteria.append("sister_clades")
                    
                    print(f"  Found cluster: {species} with {len(seq_ids)} sequences ({', '.join(criteria)})")
    
    return euk_clusters

def filter_sequences(tree, clusters, seq_lengths, output_remove_file):
    """
    Filter sequences: keep longest from each cluster, mark others for removal.
    For species with multiple clusters, check if they should be merged based on 
    phylogenetic distance.
    """
    to_remove = set()
    to_keep = set()
    
    # First, handle species with single clusters normally
    # Then handle species with multiple clusters specially
    
    for species, cluster_list in clusters.items():
        print(f"\nSpecies: {species}")
        print(f"  Found {len(cluster_list)} clusters")
        
        if len(cluster_list) == 1:
            # Single cluster - process normally
            clade, seq_ids = cluster_list[0]
            if len(seq_ids) < 2:
                to_keep.update(seq_ids)
                continue
                
            print(f"  Cluster 1: {len(seq_ids)} sequences")
            
            # Find longest sequence in this cluster
            longest_id = None
            longest_len = 0
            
            for seq_id in seq_ids:
                if seq_id not in seq_lengths:
                    print(f"    WARNING: {seq_id} not found in fasta file, skipping")
                    continue
                seq_len = seq_lengths[seq_id]
                if seq_len > longest_len:
                    longest_len = seq_len
                    longest_id = seq_id
            
            if longest_id:
                to_keep.add(longest_id)
                print(f"    Keeping: {longest_id} (length: {longest_len})")
                
                # Mark others for removal
                for seq_id in seq_ids:
                    if seq_id != longest_id and seq_id in seq_lengths:
                        to_remove.add(seq_id)
                        print(f"    Removing: {seq_id} (length: {seq_lengths[seq_id]})")
        
        else:
            # Multiple clusters - check if they should be merged
            all_seq_ids = []
            all_clades = []
            for clade, seq_ids in cluster_list:
                all_seq_ids.extend(seq_ids)
                all_clades.append(clade)
            
            print(f"  Multiple clusters with {len(all_seq_ids)} total sequences")
            
            # Find the common ancestor of all clusters for this species
            try:
                # Get terminal nodes for all sequences
                target_terminals = []
                for terminal in tree.get_terminals():
                    if terminal.name in all_seq_ids:
                        target_terminals.append(terminal)
                
                if len(target_terminals) >= 2:
                    # Find common ancestor
                    mrca = tree.common_ancestor(target_terminals)
                    mrca_descendants = get_all_descendants(tree, mrca)
                    mrca_euk_descendants = [seq_id for seq_id in mrca_descendants if is_eukaryotic_id(seq_id)]
                    
                    # Check if this species dominates its common ancestor clade
                    species_fraction = len(all_seq_ids) / len(mrca_euk_descendants)
                    
                    print(f"    Common ancestor contains {len(mrca_euk_descendants)} euk sequences")
                    print(f"    This species fraction: {species_fraction:.2f}")
                    
                    # If species dominates (>40%) or it's a small clade (≤8 euk sequences),
                    # treat all clusters as one big cluster
                    if species_fraction > 0.4 or len(mrca_euk_descendants) <= 8:
                        print(f"    Merging all clusters - species dominates common ancestor")
                        
                        # Find longest sequence across all clusters
                        longest_id = None
                        longest_len = 0
                        
                        for seq_id in all_seq_ids:
                            if seq_id not in seq_lengths:
                                print(f"      WARNING: {seq_id} not found in fasta file, skipping")
                                continue
                            seq_len = seq_lengths[seq_id]
                            if seq_len > longest_len:
                                longest_len = seq_len
                                longest_id = seq_id
                        
                        if longest_id:
                            to_keep.add(longest_id)
                            print(f"    Keeping: {longest_id} (length: {longest_len}) from merged clusters")
                            
                            # Remove all others
                            for seq_id in all_seq_ids:
                                if seq_id != longest_id and seq_id in seq_lengths:
                                    to_remove.add(seq_id)
                                    print(f"    Removing: {seq_id} (length: {seq_lengths[seq_id]})")
                    else:
                        print(f"    Keeping clusters separate - species does not dominate common ancestor")
                        # Process each cluster independently
                        for i, (clade, seq_ids) in enumerate(cluster_list):
                            print(f"  Cluster {i+1}: {len(seq_ids)} sequences")
                            
                            longest_id = None
                            longest_len = 0
                            
                            for seq_id in seq_ids:
                                if seq_id not in seq_lengths:
                                    continue
                                seq_len = seq_lengths[seq_id]
                                if seq_len > longest_len:
                                    longest_len = seq_len
                                    longest_id = seq_id
                            
                            if longest_id:
                                to_keep.add(longest_id)
                                print(f"    Keeping: {longest_id} (length: {longest_len})")
                                
                                for seq_id in seq_ids:
                                    if seq_id != longest_id and seq_id in seq_lengths:
                                        to_remove.add(seq_id)
                                        print(f"    Removing: {seq_id} (length: {seq_lengths[seq_id]})")
                
            except Exception as e:
                print(f"    Error analyzing common ancestor: {e}")
                print(f"    Processing clusters independently")
                # Fallback: process each cluster independently
                for i, (clade, seq_ids) in enumerate(cluster_list):
                    print(f"  Cluster {i+1}: {len(seq_ids)} sequences")
                    
                    longest_id = None
                    longest_len = 0
                    
                    for seq_id in seq_ids:
                        if seq_id not in seq_lengths:
                            continue
                        seq_len = seq_lengths[seq_id]
                        if seq_len > longest_len:
                            longest_len = seq_len
                            longest_id = seq_id
                    
                    if longest_id:
                        to_keep.add(longest_id)
                        print(f"    Keeping: {longest_id} (length: {longest_len})")
                        
                        for seq_id in seq_ids:
                            if seq_id != longest_id and seq_id in seq_lengths:
                                to_remove.add(seq_id)
                                print(f"    Removing: {seq_id} (length: {seq_lengths[seq_id]})")
    
    # Write only remove file
    with open(output_remove_file, 'w') as f:
        for seq_id in sorted(to_remove):
            f.write(f"{seq_id}\n")
    
    print(f"\n\nSummary:")
    print(f"  Total sequences to keep: {len(to_keep)}")
    print(f"  Total sequences to remove: {len(to_remove)}")
    print(f"  Remove list written to: {output_remove_file}")
    
    return to_keep, to_remove

def main():
    if len(sys.argv) != 4:
        print("Usage: python filter_euk_clusters.py <tree.nexus> <sequences.fasta> <remove_ids.txt>")
        print("\nThis script:")
        print("  1. Reads a phylogenetic tree in Nexus format")
        print("  2. Identifies clusters of eukaryotic sequences from the same species")
        print("  3. Keeps only the longest sequence from each cluster")
        print("  4. Outputs list of sequence IDs to remove")
        sys.exit(1)
    
    tree_file = sys.argv[1]
    fasta_file = sys.argv[2]
    remove_file = sys.argv[3]
    
    print(f"Reading tree from: {tree_file}")
    tree = Phylo.read(tree_file, "newick")
    
    print(f"Reading sequences from: {fasta_file}")
    seq_lengths = get_sequence_lengths(fasta_file)
    print(f"  Found {len(seq_lengths)} sequences")
    
    print("\nFinding eukaryotic clusters...")
    clusters = find_euk_clusters(tree)
    
    if not clusters:
        print("No eukaryotic clusters found!")
        # Create empty output file
        open(remove_file, 'w').close()
        return
    
    print(f"Found {len(clusters)} species with clustered sequences")
    
    print("\nFiltering sequences...")
    filter_sequences(tree, clusters, seq_lengths, remove_file)

if __name__ == "__main__":
    main()
