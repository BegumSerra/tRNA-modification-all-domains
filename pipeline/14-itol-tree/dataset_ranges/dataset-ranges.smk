# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

import os

configfile: "config.yaml"

PROTEIN = config["PROTEIN"]  # List of proteins
COMBINED_PROTEIN = config.get("COMBINED_PROTEIN", [])
TREE = config["TREE"]
protein_to_tree = config["MAP"]  # Dictionary mapping proteins to trees
tree_to_protein = {tree_val: protein_key for protein_key, tree_val in protein_to_tree.items()}

combined_fasta_base = config.get("fasta_combined_tree")
if combined_fasta_base and not os.path.exists(combined_fasta_base):
    combined_fasta_base = combined_fasta_base.replace("combined-tree", "combined_tree")
if not combined_fasta_base:
    combined_fasta_base = config["fasta"].replace("homologs", "combined_tree")

wildcard_constraints:
    protein = "|".join(PROTEIN),
    combined = "|".join(COMBINED_PROTEIN)

def get_protein_id_from_tree(tree):
    if tree not in tree_to_protein:
        raise ValueError(f"tree '{tree}' not found in mapping.")
    return tree_to_protein[tree]

def get_tree(wildcards):
    return f"{config['tree']}/{protein_to_tree[wildcards.protein]}_Cmod.treefile"

def get_fasta(wildcards):
    return f"{config['fasta']}/{protein_to_tree[wildcards.protein]}.fa"

def get_tree_combined(wildcards):
    return f"{config['tree']}/combined_tree/{protein_to_tree[wildcards.combined]}_Cmod.treefile"

def get_fasta_combined(wildcards):
    return f"{combined_fasta_base}/{protein_to_tree[wildcards.combined]}.fa"

def get_protein_arg(wildcards):
    return wildcards.protein

rule all:
    input:
         expand(config["label"] +"/{protein}_labels.txt", protein=PROTEIN),
         expand(config["label"] +"/combined_tree/{combined}_labels.txt", combined=COMBINED_PROTEIN),
         expand(config["range"] +"/{protein}_choloro_euk.txt", protein=[p for p in PROTEIN if p in ["MnmG", "MnmE", "TilS"]]),
         expand(config["range"] +"/{protein}_dataset_ranges.txt", protein=PROTEIN),
         expand(config["range"] +"/combined_tree/{combined}_dataset_ranges.txt", combined=COMBINED_PROTEIN),
         expand(config["range"] +"/combined_tree/{proteins}_protein_layer.txt", proteins=COMBINED_PROTEIN),
         expand(config["range"] +"/combined_tree/{proteins}_protein_text.txt", proteins=COMBINED_PROTEIN)

rule create_labels:
    input:
        fasta = get_fasta
    output:
        label = config["label"] +"/{protein}_labels.txt"
    shell:
        "python fasta_to_csv.py {input.fasta} {output.label}"

rule create_labels_combined:
    input:
        fasta = get_fasta_combined
    output:
        label = config["label"] +"/combined_tree/{combined}_labels.txt"
    shell:
        "python fasta_to_csv.py {input.fasta} {output.label}"

rule find_mrca:
    input:
        tree = get_tree
    output:
        ranges = config["range"] +"/{protein}_choloro_euk.txt"
    params:
        protein = get_protein_arg
    wildcard_constraints:
        protein = "MnmG|MnmE|TilS"
    shell:
        "python find_mrca.py {input.tree} {output.ranges} {params.protein}"

rule dataset_ranges_with_mrca:
    input:
        fasta = get_fasta,
        phylum = config["phylum_table"] + "/{protein}_phylum.table.tsv",
        mrca_ids = config["range"] + "/{protein}_choloro_euk.txt"
    output:
        ranges = config["range"] +"/{protein}_dataset_ranges.txt"
    params:
        protein = get_protein_arg
    wildcard_constraints:
        protein = "MnmG|MnmE|TilS"
    shell:
        "python dataset-ranges.py {input.fasta} {input.phylum} {output.ranges} {params.protein} {input.mrca_ids}"

rule dataset_ranges:
    input:
        fasta = get_fasta,
        phylum = config["phylum_table"] + "/{protein}_phylum.table.tsv"
    output:
        ranges = config["range"] +"/{protein}_dataset_ranges.txt"
    params:
        protein = get_protein_arg
    wildcard_constraints:
        protein = "MnmH|TgtA|Tuc1|TiaS|ArcS|Elp3"
    shell:
        "python dataset-ranges.py {input.fasta} {input.phylum} {output.ranges} {params.protein}"

rule dataset_ranges_combined:
    input:
        fasta = get_fasta_combined,
        phylum = config["phylum_table"] + "/combined_tree/{combined}_phylum.table.tsv"
    output:
        ranges = config["range"] +"/combined_tree/{combined}_dataset_ranges.txt"
    params:
        protein = lambda wildcards: wildcards.combined
    wildcard_constraints:
        combined = "|".join(COMBINED_PROTEIN) if COMBINED_PROTEIN else "^$"
    shell:
        "python dataset-ranges.py {input.fasta} {input.phylum} {output.ranges} {params.protein}"

# Mapping for individual proteins that are part of combined trees
# These map protein names to their FASTA file basenames (without .fa)
# NOTE: Use pre-r01 files (nogaps without _r01) because combined trees are built from 
# pre-r01 files and the r01 filtering happens AFTER merging. Using post-r01 files would 
# cause sequences in the combined tree to be missing from the protein layer annotation.
INDIVIDUAL_PROTEIN_MAP = config.get("INDIVIDUAL_PROTEIN_MAP", {
    "MnmA": "MnmA_k01_f01_e01_c01_nogaps",
    "MTU1": "MTU1_k01_f01_e00_c01_nogaps",
    "TadA": "TadA_k01_f01_e01_c01_nogaps",
    "ADAT2": "ADAT2_k01_f01_e01_c01_nogaps",
    "Tgt": "Tgt_k01_f00_e01_c01_nogaps",
    "QTRT1": "QTRT1_k01_f00_e00_c01_nogaps",
})

def get_individual_fastas_for_combined(wildcards):
    """Get individual FASTA files for each protein in a combined tree (e.g., TadA-ADAT2 -> TadA.fa, ADAT2.fa)"""
    proteins = wildcards.combined.split("-")
    fasta_files = []
    for p in proteins:
        if p in INDIVIDUAL_PROTEIN_MAP:
            fasta_files.append(f"{config['fasta']}/{INDIVIDUAL_PROTEIN_MAP[p]}.fa")
        elif p in protein_to_tree:
            fasta_files.append(f"{config['fasta']}/{protein_to_tree[p]}.fa")
        else:
            raise ValueError(f"Protein '{p}' not found in INDIVIDUAL_PROTEIN_MAP or protein_to_tree")
    return fasta_files

rule dataset_ranges_protein_layer:
    input:
        fasta_files = get_individual_fastas_for_combined
    output:
        ranges = config["range"] +"/combined_tree/{combined}_protein_layer.txt"
    wildcard_constraints:
        combined = "|".join(COMBINED_PROTEIN) if COMBINED_PROTEIN else "^$"
    shell:
        "python dataset-ranges-protein-layer.py {input.fasta_files} {output.ranges}"

rule dataset_text_protein_layer:
    input:
        fasta_files = get_individual_fastas_for_combined
    output:
        text = config["range"] +"/combined_tree/{combined}_protein_text.txt"
    wildcard_constraints:
        combined = "|".join(COMBINED_PROTEIN) if COMBINED_PROTEIN else "^$"
    shell:
        "python dataset-text-protein-layer.py {input.fasta_files} {output.text}"
