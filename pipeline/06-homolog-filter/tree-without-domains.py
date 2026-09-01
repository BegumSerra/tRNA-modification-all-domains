# Copyright (c) 2024, Begüm Serra Büyüktarakçı
# All rights reserved.

import sys
import pandas as pd
from Bio import SeqIO
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, AttrFace


def load_protein_mapping(mapping_file):
    """Load protein-genome mapping to get species names and domains"""
    df = pd.read_csv(mapping_file, sep='\t')
    # Create a dictionary: protein_id -> species_name and domain
    # Use drop_duplicates to keep first occurrence if protein appears in multiple genomes
    df_unique = df.drop_duplicates(subset='protein_id', keep='first')
    mapping = {}
    for _, row in df_unique.iterrows():
        mapping[row['protein_id']] = {
            'species': row['species'],
            'domain': row['domain']
        }
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
    tree_file, merged_fasta_file, mapping_file, output_image = sys.argv[1:5]

    # Load protein-genome mapping
    protein_mapping = load_protein_mapping(mapping_file)

    # Process the tree for annotations
    t = process_tree(tree_file, protein_mapping)

    # Ladderize the tree
    t.ladderize()
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)

    # Set display options and render the tree
    ts = display_options()
    t.render(output_image, tree_style=ts)
