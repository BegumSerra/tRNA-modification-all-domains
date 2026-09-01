# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "tree.yaml"
results = "../../output/"
data = "../../input/"
envs = "../../envs/"
database = "../../input/database/"
PROTEIN = config["PROTEIN"]  # List of filter strings
filter_to_protein = config["FILTER_TO_PROTEIN"]  # protein_name -> filter_string
# Create reverse mapping: filter_string -> protein_name
protein_to_filter = {v: k for k, v in filter_to_protein.items()}

# For wildcards, we use protein names (extracted from filter strings)
PROTEIN = [p for p in filter_to_protein.keys() if p != "TusC"]

wildcard_constraints:
    protein="|".join(PROTEIN)

rule all:
    input:
        expand(results + "tree_images/{protein}_nogaps.domain.tree.pdf", protein=PROTEIN),
        expand(results + "tree_images/{protein}_nogaps.no_domain.tree.pdf", protein=PROTEIN)

rule with_domain:
    input:
        domain_file = results + "domain/blast/merged/{protein}.tsv",
        phylum_table = results + "table/phylum_table/{protein}_phylum.table.tsv",
        tree_file = lambda wildcards: results + f"iqtree/{filter_to_protein[wildcards.protein]}_nogaps.treefile",
        fasta_file = lambda wildcards: results + f"seqs/homologs/{filter_to_protein[wildcards.protein]}_nogaps.fa"
    output:
        results + "tree_images/{protein}_nogaps.domain.tree.pdf"
    shell:
        "conda run -n ete python tree-with-domains.py {input.domain_file} {input.phylum_table} {input.tree_file} {input.fasta_file} {output}"

rule without_domain:
    input:
        tree = lambda wildcards: results + f"iqtree/{filter_to_protein[wildcards.protein]}_nogaps.treefile",
        fasta = lambda wildcards: results + f"seqs/homologs/{filter_to_protein[wildcards.protein]}_nogaps.fa",
        mapping = results + "table/phylum_table/{protein}_phylum.table.tsv"
    output:
        results + "tree_images/{protein}_nogaps.no_domain.tree.pdf"
    shell:
        "conda run -n ete python tree-without-domains.py {input.tree} {input.fasta} {input.mapping} {output}"
