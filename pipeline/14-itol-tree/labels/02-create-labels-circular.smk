# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "../../../config/28-itol.yaml"
wildcard_constraints:
    protein="|".join(config["PROTEIN"])
results = "../../../output/"
data = "../../../input/"
envs = "../../../envs/"
database = "../../../database/"
# Load the protein and tree mapping from the config file
PROTEIN = config["PROTEIN"]  # List of proteins
TREE = config["TREE"]
protein_to_tree = config["MAP"]  = config["MAP"]  # Dictionary mapping proteins to trees
tree_to_protein = {tree_val: protein_key for protein_key, tree_val in protein_to_tree.items()}

def get_protein_id_from_tree(tree):
    if tree not in tree_to_protein:
        raise ValueError(f"Tree '{tree}' not found in mapping.")
    return tree_to_protein[tree]

rule all:
    input:
        expand(results + "itol/templates/labels/data_{protein}.txt", protein=PROTEIN),
        expand(results + "itol/templates/labels/labels_template_{protein}.txt", protein=PROTEIN)

rule create_labels:
    input:
        lambda wildcards: results + f"seqs/filtering/{protein_to_tree[wildcards.protein]}.fa"
    output:
        results + "itol/templates/labels/data_{protein}.txt"
    shell:
        "python short-species-name-labels.py {input} {output}"
