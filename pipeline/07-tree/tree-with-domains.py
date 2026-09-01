# coding=utf-8
# Got from Anna Odelgard and modified by Begüm Serra Büyüktarakçı 2021

## Run on command line:
## python domain_tree_clade_flag.py /Users/serbu540/Projects/project_first/results/domain/RluABD_merged.fa.tsv /Users/serbu540/Projects/project_first/results/phylum_table/RluABD_phylum_table.tsv /Users/serbu540/Projects/project_first/data/tables/genome_statistics.tsv /Users/serbu540/Projects/project_first/results/align_sco/iqtree_files/RluABD.contree /Users/serbu540/Projects/project_first/results/seqs/RluABD_merged.fa /Users/serbu540/Projects/project_first/results/align_sco/tree_images/RLUABD_domain.png

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
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace, NodeStyle
from Bio import SeqIO
import argparse

def extract_domains(domain_file):
    domain_df = domain_file[domain_file["Analysis"] == "Pfam"]
    return domain_df

def retrieve_pids(fasta_file_path):
    """
    retrieve the protein ids from the fasta file
    """
    fasta_file = SeqIO.parse(fasta_file_path, "fasta")
    protein_ids = []

    for s in fasta_file:
        protein_ids.append(s.id)

    return protein_ids

def domain_info(protein_ids, domain_df):
    """
    create a list of domain motifs to be added to the tree
    """
    domain_info_list = []
    domain_info_list_temp = []
    color = set_color(domain_df)

    # Go trough the protein ids and extract the domain records in the intersproscan file
    for pid in protein_ids:
        protein_domains = domain_df[domain_df['Protein accession']==pid]
        protein_domains = protein_domains.reset_index(drop=True)
        #print(protein_domains)
        # for each record in the interproscan file add the domain info to the domain info list
        for row in protein_domains.index:

            domain = protein_domains.iloc[row,:]
            domain_start = domain['Start location']
            domain_stop = domain['Stop location']
            domain_id = domain['Signature accession']
            motif = [domain_start, domain_stop, "()", None, 14, "black", color[domain_id], "arial|10|black|" + str(domain_id)]
            domain_info_list_temp.append(motif)
            motif = []
        domain_info_list.append([pid,domain_info_list_temp])
        domain_info_list_temp = []
    #print(domain_info_list)

    return domain_info_list

def set_color(domain_df):

    color = {}
    color_index = 0
    color_list = [
        '#01cdfe', '#ff71ce', '#05ffa1', '#b967ff', '#fffb96', '#ffa105', '#ffacfb', '#5e9bff',
        '#ff5f5f', '#09A67B', '#d5a3ff', '#954d5d', '#b74f6f', '#cd6c6c', '#dd8a6e', '#e6aa68',
        '#4bdebf', '#f24830', '#fefe14', '#ff98c7', '#0eabce', '#ff69b4', '#8c9eff', '#ffc300',
        '#7fffd4', '#ff6347', '#40e0d0', '#fa8072', '#4682b4', '#9370db',
        '#3cb44b', '#ffd700', '#ff4500', '#6495ed', '#dc143c', '#8b4513', '#00ced1',
        '#ffdab9', '#ee82ee', '#bc8f8f',
        '#01cdfe', '#ff71ce', '#05ffa1', '#b967ff', '#fffb96', '#ffa105', '#ffacfb', '#5e9bff',
        '#ff5f5f', '#09A67B', '#d5a3ff', '#954d5d', '#b74f6f', '#cd6c6c', '#dd8a6e', '#e6aa68',
        '#4bdebf', '#f24830', '#fefe14', '#ff98c7', '#0eabce', '#ff69b4', '#8c9eff', '#ffc300',
        '#7fffd4', '#ff6347', '#40e0d0', '#fa8072', '#4682b4', '#9370db',
        '#3cb44b', '#ffd700', '#ff4500', '#6495ed', '#dc143c', '#8b4513', '#00ced1',
        '#ffdab9', '#ee82ee', '#bc8f8f',
        '#01cdfe', '#ff71ce', '#05ffa1', '#b967ff', '#fffb96', '#ffa105', '#ffacfb', '#5e9bff',
        '#ff5f5f', '#09A67B', '#d5a3ff', '#954d5d', '#b74f6f', '#cd6c6c', '#dd8a6e', '#e6aa68',
        '#4bdebf', '#f24830', '#fefe14', '#ff98c7', '#0eabce', '#ff69b4', '#8c9eff', '#ffc300',
        '#7fffd4', '#ff6347', '#40e0d0', '#fa8072', '#4682b4', '#9370db',
        '#3cb44b', '#ffd700', '#ff4500', '#6495ed', '#dc143c', '#8b4513', '#00ced1',
        '#ffdab9', '#ee82ee', '#bc8f8f',
        '#01cdfe', '#ff71ce', '#05ffa1', '#b967ff', '#fffb96', '#ffa105', '#ffacfb', '#5e9bff',
        '#ff5f5f', '#09A67B', '#d5a3ff', '#954d5d', '#b74f6f', '#cd6c6c', '#dd8a6e', '#e6aa68',
        '#4bdebf', '#f24830', '#fefe14', '#ff98c7', '#0eabce', '#ff69b4', '#8c9eff', '#ffc300',
        '#7fffd4', '#ff6347', '#40e0d0', '#fa8072', '#4682b4', '#9370db',
        '#3cb44b', '#ffd700', '#ff4500', '#6495ed', '#dc143c', '#8b4513', '#00ced1',
        '#ffdab9', '#ee82ee', '#bc8f8f'
    ]

    domain_ids_list = domain_df['Signature accession'].drop_duplicates().values.tolist()
    #print(domain_ids_list)
    for domain_id in domain_ids_list:
        color[domain_id]= color_list[color_index]

        #print(color[domain_id])
        color_index+=1

    return color

# add changes to the tree, add domain figure and organism names
def add_tree_info(t, fasta_file_path, domain_info_list):
    for node in t.traverse():
        sequence = []
        motif = []
        if node.is_leaf():
            accession = node.name

            # retrieve the fasta sequence for the protein
            fasta_file = SeqIO.parse(fasta_file_path, "fasta")

            for s in fasta_file:
                if s.id.rstrip() == accession.rstrip():
                    sequence = s.seq

            # retrieve the motif for the protein
            for line in domain_info_list:
                if line[0] == node.name:
                    #print(node.name)
                    motif = line[1]
                    #print(motif)

            # add the information to the tree
            SeqFace = SeqMotifFace(str(sequence), motifs = motif, seq_format="-")
            (t & node.name).add_face(SeqFace, 0, "aligned")

    return t

def style_tree(tree, protein_mapping):
    """
    Add color and organism name to taxa
    """

    colours = {'Halobacteriota': '#CC79A7',
              'Methanobacteriota' : '#70000E',
              'Thermoplasmatota': '#00702a',
              'Methanobacteriota_B': '#ec5f67',
              'Thermoproteota' : '#6699cc',
              'Archaea': '#648FFF',
              'Bacteria': '#DC267F',
              'Eukarya': '#FE6100'}

    nstyle_supported = NodeStyle()
    nstyle_supported['size'] = 5
    nstyle_supported['fgcolor'] = 'black'
    nstyle_unsupported = NodeStyle()
    nstyle_unsupported['size'] = 0

    for node in tree.traverse():
            # if node support >=90 put black dot
            if node.support >= 90 and node.is_leaf() == False:
                node.set_style(nstyle_supported)
            else:
                node.set_style(nstyle_unsupported)

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
                
                # Check if protein exists in mapping
                if protein_id not in protein_mapping:
                    print(f"WARNING: Protein {protein_id} not found in mapping, skipping styling")
                    continue
                
                species = protein_mapping[protein_id]['species']
                domain = protein_mapping[protein_id]['domain']
                print(protein_id)

                #change the acc number on the tree to species name and color
                organism_name = TextFace(species, fgcolor=colours.get(domain, 'black'))
                organism_name.margin_left = 10
                node.add_face(organism_name, column=0, position='branch-right')
    return tree

def ladarize_tree(tree):
    tree.ladderize()
    R = tree.get_midpoint_outgroup()
    # and set it as tree outgroup
    tree.set_outgroup(R)
    return tree

def display_options(tree):
    ts = TreeStyle()
    ts.show_branch_support = True
    ts.show_branch_length = False
    ts.complete_branch_lines_when_necessary = True  # Connect branches without gaps
    ts.show_leaf_name = True
    ts.scale = 100
    #ts.title.add_face(TextFace(figure_name, fsize=20), column=0)
    return ts

###flags###
parser = argparse.ArgumentParser()
parser.add_argument("in_domainfile", type=str, help="domain file path")
parser.add_argument("in_phylumtable", type=str, help="phylum table path")
parser.add_argument("in_raw_treefile", type=str, help="raw tree file path")
parser.add_argument("in_fastafile", type=str, help="fasta file path")
#parser.add_argument("in_OG_table", type=str, help="OG table path")

parser.add_argument("out_edited_treefile", type=str, help="fasta file with annotation")
args = parser.parse_args()

#input
domain_file_path = args.in_domainfile
mapping_file = args.in_phylumtable  # Now points to protein_genome_mapping.tsv
tree_path = args.in_raw_treefile
fasta_file_path = args.in_fastafile
#OG_table = args.in_OG_table

#my_OG = str(tree_path).split("/")[-1].split("_")[0]
#output
file_name = args.out_edited_treefile

#domain_file
interpro_header = ['Protein accession','Sequence MD5 digest','Sequence length','Analysis','Signature accession','Signature description','Start location','Stop location','Score','Status','Date','InterPro accession','InterPro description','GO annotations','Pathways annotations']
domain_file= pd.read_csv(args.in_domainfile, sep='\t', names = interpro_header)

# Load protein-genome mapping
mapping_df = pd.read_csv(mapping_file, sep='\t')
# Drop duplicates to handle proteins in multiple genomes
mapping_df = mapping_df.drop_duplicates(subset='acc', keep='first')
# Create dictionary for fast lookup
protein_mapping = {}
for _, row in mapping_df.iterrows():
    protein_mapping[row['acc']] = {
        'species': row['species'],
        'domain': row['domain']
    }

# Load tree
tree = Tree(tree_path)

#to run functions:
domain_df = extract_domains(domain_file)
protein_ids = retrieve_pids(fasta_file_path)
domain_info_list = domain_info(protein_ids, domain_df)
color= set_color(domain_df)
t = style_tree(tree, protein_mapping)
t = ladarize_tree(t)
ts = display_options(t)
t = add_tree_info(t, fasta_file_path, domain_info_list)
try:
    t.render(file_name, tree_style=ts, w=4000)
    print(f"Tree rendered successfully to {file_name}")
except Exception as e:
    print(f"Render error: {e}")
import sys
sys.exit(0)
#tree.show()
