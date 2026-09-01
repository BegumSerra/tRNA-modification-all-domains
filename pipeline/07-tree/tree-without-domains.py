# Copyright (c) 2024, Begüm Serra Büyüktarakçı
# All rights reserved.

import sys
import os

# CRITICAL: Set these BEFORE any other imports to prevent Qt initialization
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
os.environ['MPLBACKEND'] = 'Agg'
os.environ['ETS_TOOLKIT'] = 'null'

# Force matplotlib to use non-interactive backend before importing ete3
import matplotlib
matplotlib.use('Agg')

import pandas as pd
from Bio import SeqIO

# Disable Qt for ete3 by importing it last
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, AttrFace


def load_protein_mapping(mapping_file):
    """Load phylum table file with columns: acc, domain, phylum, order, family, species"""
    print(f"Loading mapping from {mapping_file}...")
    df = pd.read_csv(mapping_file, sep='\t')
    
    mapping = {}
    for _, row in df.iterrows():
        # acc column contains protein_id
        protein_id = row['acc'] if pd.notna(row['acc']) else ''
        if protein_id:  # Skip empty rows
            mapping[protein_id] = {
                'species': row['species'] if pd.notna(row['species']) else '',
                'domain': row['domain'] if pd.notna(row['domain']) else ''
            }
    
    print(f"Loaded {len(mapping)} proteins")
    return mapping

def fasta_to_dataframe(fasta_file):
    data = {'id': [], 'description': []}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        data['id'].append(record.id)
        data['description'].append(record.description)
    return pd.DataFrame(data)

def process_tree(tree_file, protein_mapping):
    """
    Add color and organism name to taxa based on domain
    """
    colours = {'Archaea': '#648FFF',
              'Bacteria': '#DC267F',
              'Eukarya': '#FE6100'}
    
    # Load the tree
    t = Tree(tree_file)

    # Node styles for branch support visualization
    nstyle_supported = NodeStyle()
    nstyle_supported['size'] = 5
    nstyle_supported['fgcolor'] = 'black'
    nstyle_unsupported = NodeStyle()
    nstyle_unsupported['size'] = 0

    # Traverse the tree to set styles and annotations
    for node in t.traverse():
        # if node support >=90 put black dot (only for non-leaf nodes)
        if node.support >= 90 and node.is_leaf() == False:
            node.set_style(nstyle_supported)
        else:
            node.set_style(nstyle_unsupported)

        # Set annotations for leaf nodes
        if node.is_leaf():
            protein_id = node.name
            
            # Check if it's a eukaryotic ID format (EP#####_Genus_species_P######)
            if protein_id.startswith('EP') and '_' in protein_id:
                parts = protein_id.split('_')
                if len(parts) >= 3:
                    # Extract genus and species from ID
                    genus = parts[1] if len(parts) > 1 else ''
                    species_epithet = parts[2] if len(parts) > 2 else ''
                    species = f"{genus} {species_epithet}".replace('_', ' ')
                    domain = 'Eukarya'
                    
                    # Add species name with color for Eukarya
                    organism_name = TextFace(species, fgcolor=colours.get(domain, 'black'))
                    organism_name.margin_left = 10
                    node.add_face(organism_name, column=0, position='branch-right')
                    print(f"{protein_id} -> {species} (Eukarya)")
                    continue
            
            if protein_id in protein_mapping:
                species = protein_mapping[protein_id]['species']
                domain = protein_mapping[protein_id]['domain']
                print(protein_id)
                
                # Add species name with color based on domain
                organism_name = TextFace(species, fgcolor=colours.get(domain, 'black'))
                organism_name.margin_left = 10
                node.add_face(organism_name, column=0, position='branch-right')
            else:
                print(f"WARNING: Protein {protein_id} not found in mapping")

    return t

def ladderize_tree(t):
    t.ladderize()
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    return t


def display_options():
    ts = TreeStyle()
    ts.show_branch_support = True
    ts.show_branch_length = False
    ts.complete_branch_lines_when_necessary = False  # Connect branches without gaps
    ts.show_leaf_name = True
    ts.scale = 100
    return ts

if __name__ == '__main__':
    try:
        tree_file, merged_fasta_file, mapping_file, output_image = sys.argv[1:5]
        
        print(f"Processing tree: {tree_file}")
        print(f"Output: {output_image}")

        # Load phylum table file
        protein_mapping = load_protein_mapping(mapping_file)

        # Process the tree for annotations
        print("Processing tree...")
        t = process_tree(tree_file, protein_mapping)

        # Ladderize the tree
        print("Ladderizing tree...")
        t.ladderize()
        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)

        # Set display options and render the tree
        print("Rendering tree...")
        ts = display_options()
        t.render(output_image, tree_style=ts)
        
        print(f"Successfully created: {output_image}")
    
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    t = process_tree(tree_file, protein_mapping)

    # Ladderize the tree
    t.ladderize()
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)

    # Set display options and render the tree
    ts = display_options()
    t.render(output_image, tree_style=ts)
