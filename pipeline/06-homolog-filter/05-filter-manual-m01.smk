# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "05-config.yaml"
wildcard_constraints:
    protein = "|".join(config["PROTEIN"]),
    filter = "|".join(config["FILTER"])

results = "../../output/"
data = "../../input/"
database = "../../database/"
# Load the protein and tree mapping from the config file
PROTEIN = config["PROTEIN"]  # List of proteins
FILTER = config["FILTER"]
filter_to_protein = config["MAP"]  # Dictionary mapping filter names to protein names

def get_protein_id_from_filter(filter_name):
    if filter_name not in filter_to_protein:
        raise ValueError(f"Filter '{filter_name}' not found in mapping.")
    return filter_to_protein[filter_name]

rule all:
    input:
        expand(results + "seqs/filtering/fasta/{filter}_m01.fa",filter=FILTER),
        #expand(results + "alignments/raw_alignments/{filter}_m01.aln", filter=FILTER),
        #expand(results+ "alignments/trimmed_alignments/{filter}_m01.trimal", filter=FILTER),
        #expand(results + "iqtree/{filter}_m01.treefile", filter=FILTER),
        #expand(results + "tree_images/{filter}_m01.domain.tree.pdf", filter=FILTER)

rule remove:
    input:
        fasta = results + "seqs/filtering/fasta/{filter}.fa",
        ids = lambda wildcards: results + f"seqs/filtering/ids/m01/{get_protein_id_from_filter(wildcards.filter)}.txt"
    output:
         results + "seqs/filtering/fasta/{filter}_m01.fa"
    shell:
        "python fasta-remove.py {input.fasta} {input.ids} {output}"

rule msa:
    input:
        results + "seqs/filtering/fasta/{filter}_m01.fa"
    output:
        results + "alignments/raw_alignments/{filter}_m01.aln"
    shell:
        "mafft --auto {input} > {output}"

rule trim:
    input:
        results + "alignments/raw_alignments/{filter}_m01.aln"
    output:
        results + "alignments/trimmed_alignments/{filter}_m01.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule tree:
    """
    generate the phylogenies using iqtree
    """
    input:
        results + "alignments/trimmed_alignments/{filter}_m01.trimal"
    output:
        results + "iqtree/{filter}_m01.treefile"
    params:
        uf=config["uf"],
        filename= results + "iqtree/{filter}_m01"
    shell:
        #"iqtree -s {input} -m MFP -mset LG -bb {params.uf} -nt {threads} -ntmax {threads} --prefix {params.filename} -redo"
        "iqtree -s {input} -m MFP -mset LG -bb {params.uf} --prefix {params.filename} -redo"

rule with_domain:
    input:
        lambda wildcards: results + f"domain/blast/merged/{get_protein_id_from_filter(wildcards.filter)}.tsv",
        lambda wildcards: results + f"table/phylum_table/{get_protein_id_from_filter(wildcards.filter)}_phylum.table.tsv",
        results + "iqtree/{filter}_m01.treefile",
        results + "seqs/filtering/fasta/{filter}_m01.fa"
    output:
        results + "tree_images/{filter}_m01.domain.tree.pdf"
    shell:
        "python tree-with-domains.py {input} {output}"
