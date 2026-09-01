# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "01-config.yaml"
results = "../../../output/"
data = "../../../input/"
database = "../../../database/"
# Load the protein and tree mapping from the config file
PROTEIN = config["PROTEIN"]  # List of regular proteins
COMBINED_PROTEIN = config.get("COMBINED_PROTEIN", [])  # List of combined proteins
TREE = config["TREE"]
COMBINED_TREE = config.get("COMBINED_TREE", [])
protein_to_tree = config["MAP"]  # Dictionary mapping proteins to trees
tree_to_protein = {tree_val: protein_key for protein_key, tree_val in protein_to_tree.items()}

wildcard_constraints:
    protein="|".join(PROTEIN),
    combined="|".join(COMBINED_PROTEIN)

def get_protein_id_from_tree(tree):
    if tree not in tree_to_protein:
        raise ValueError(f"Tree '{tree}' not found in mapping.")
    return tree_to_protein[tree]

rule all:
    input:
        expand(results + "itol/templates/labels/data_{protein}.txt", protein=PROTEIN),
        expand(results + "itol/templates/labels/text_{protein}.txt", protein=PROTEIN),
        expand(results + "itol/templates/labels/combined_tree/data_{combined}.txt", combined=COMBINED_PROTEIN),
        expand(results + "itol/templates/labels/combined_tree/text_{combined}.txt", combined=COMBINED_PROTEIN)

rule create_labels:
    input:
        fasta = lambda wildcards: results + f"seqs/homologs/{protein_to_tree[wildcards.protein]}.fa",
        mapping = results + "table/protein_genome_mapping.tsv"
    output:
        results + "itol/templates/labels/data_{protein}.txt"
    shell:
        "python ids-long-species-name-labels.py {input.fasta} {input.mapping} {output}"

rule create_text_labels:
    input:
        fasta = lambda wildcards: results + f"seqs/homologs/{protein_to_tree[wildcards.protein]}.fa",
        mapping = results + "table/protein_genome_mapping.tsv"
    output:
        results + "itol/templates/labels/text_{protein}.txt"
    shell:
        "python create-text-labels.py {input.fasta} {input.mapping} {output}"

# Rules for combined tree proteins
rule create_labels_combined:
    input:
        fasta = lambda wildcards: results + f"seqs/combined_tree/{protein_to_tree[wildcards.combined]}.fa",
        mapping = results + "table/protein_genome_mapping.tsv"
    output:
        results + "itol/templates/labels/combined_tree/data_{combined}.txt"
    shell:
        "python ids-long-species-name-labels.py {input.fasta} {input.mapping} {output}"

rule create_text_labels_combined:
    input:
        fasta = lambda wildcards: results + f"seqs/combined_tree/{protein_to_tree[wildcards.combined]}.fa",
        mapping = results + "table/protein_genome_mapping.tsv"
    output:
        results + "itol/templates/labels/combined_tree/text_{combined}.txt"
    shell:
        "python create-text-labels.py {input.fasta} {input.mapping} {output}"
