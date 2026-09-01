# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "edit_all.yaml"
results = "../../output/"
data = "../../input/"
envs = "../../envs/"
database = "../../input/database/"
SINGLE = config.get("SINGLE", [])  # List of single protein filter strings
COMBINED = config.get("COMBINED", [])  # List of combined protein filter strings
filter_to_protein = config["FILTER_TO_PROTEIN"]  # protein_name -> filter_string
# Create reverse mapping: filter_string -> protein_name
protein_to_filter = {v: k for k, v in filter_to_protein.items()}

# Separate single and combined proteins from FILTER_TO_PROTEIN
SINGLE_PROTEINS = [p for p in filter_to_protein.keys() if '-' not in p and p != "TusC"]
COMBINED_PROTEINS = [p for p in filter_to_protein.keys() if '-' in p]

wildcard_constraints:
    protein="|".join(SINGLE_PROTEINS),
    combined="|".join(COMBINED_PROTEINS)

rule all:
    input:
        expand(results + "tree_images/Cmod/{protein}.domain.Cmod.tree.pdf", protein=SINGLE_PROTEINS),
        expand(results + "tree_images/Cmod/{protein}.no_domain.Cmod.tree.pdf", protein=SINGLE_PROTEINS),
        expand(results + "tree_images/Cmod/combined_tree/{combined}.domain.Cmod.tree.pdf", combined=COMBINED_PROTEINS),
        expand(results + "tree_images/Cmod/combined_tree/{combined}.no_domain.Cmod.tree.pdf", combined=COMBINED_PROTEINS)

# Rules for single proteins (like ELP3 with m01 filtering)
rule with_domain:
    input:
        domain_file = results + "domain/blast/merged/{protein}.tsv",
        phylum_table = results + "table/phylum_table/{protein}_phylum.table.tsv",
        tree_file = lambda wildcards: results + f"iqtree/{filter_to_protein[wildcards.protein]}_Cmod.treefile",
        fasta_file = lambda wildcards: results + f"seqs/homologs/{filter_to_protein[wildcards.protein]}.fa"
    output:
        results + "tree_images/Cmod/{protein}.domain.Cmod.tree.pdf"
    shell:
        "conda run -n ete python ../07-tree/tree-with-domains.py {input.domain_file} {input.phylum_table} {input.tree_file} {input.fasta_file} {output}"

rule without_domain:
    input:
        tree = lambda wildcards: results + f"iqtree/{filter_to_protein[wildcards.protein]}_Cmod.treefile",
        fasta = lambda wildcards: results + f"seqs/homologs/{filter_to_protein[wildcards.protein]}.fa",
        mapping = results + "table/phylum_table/{protein}_phylum.table.tsv"
    output:
        results + "tree_images/Cmod/{protein}.no_domain.Cmod.tree.pdf"
    shell:
        "conda run -n ete python ../07-tree/tree-without-domains.py {input.tree} {input.fasta} {input.mapping} {output}"

# Rules for combined tree proteins
rule with_domain_combined:
    input:
        domain_file = lambda wildcards: results + f"domain/combined_tree/{wildcards.combined.split('_')[0]}.tsv",
        phylum_table = lambda wildcards: results + f"phylum_table/combined_tree/{wildcards.combined.split('_')[0]}_phylum.table.tsv",
        tree_file = lambda wildcards: results + f"iqtree/combined_tree/{filter_to_protein[wildcards.combined]}_Cmod.treefile",
        fasta_file = lambda wildcards: results + f"seqs/combined_tree/{filter_to_protein[wildcards.combined]}.fa"
    output:
        results + "tree_images/Cmod/combined_tree/{combined}.domain.Cmod.tree.pdf"
    shell:
        "conda run -n ete python ../07-tree/tree-with-domains.py {input.domain_file} {input.phylum_table} {input.tree_file} {input.fasta_file} {output}"

rule without_domain_combined:
    input:
        tree = lambda wildcards: results + f"iqtree/combined_tree/{filter_to_protein[wildcards.combined]}_Cmod.treefile",
        fasta = lambda wildcards: results + f"seqs/combined_tree/{filter_to_protein[wildcards.combined]}.fa",
        mapping = lambda wildcards: results + f"phylum_table/combined_tree/{wildcards.combined.split('_')[0]}_phylum.table.tsv"
    output:
        results + "tree_images/Cmod/combined_tree/{combined}.no_domain.Cmod.tree.pdf"
    shell:
        "conda run -n ete python ../07-tree/tree-without-domains.py {input.tree} {input.fasta} {input.mapping} {output}"

