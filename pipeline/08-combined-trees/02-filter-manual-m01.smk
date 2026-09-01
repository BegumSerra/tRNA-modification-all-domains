# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "config.yaml"

results = "../../output/"
data = "../input/"

# Load configuration
MERGE = config["MERGE"]  # List of merged protein names (e.g., Tgt-QTRT1, TadA-ADAT2)
PROTEINS_MAP = config["PROTEINS"]  # Dictionary mapping merged names to list of proteins

wildcard_constraints:
    merge = "|".join(MERGE)

rule all:
    input:
        expand(results + "seqs/combined_tree/{merge}_m01.fa", merge=MERGE),
        expand(results + "alignments/combined_tree/{merge}_m01.aln", merge=MERGE),
        expand(results + "alignments/combined_tree/{merge}_m01.trimal", merge=MERGE),
        expand(results + "iqtree/combined_tree/{merge}_m01.treefile", merge=MERGE),
        expand(results + "tree_images/combined_tree/{merge}_m01.domain.tree.pdf", merge=MERGE)

rule remove:
    input:
        fasta = results + "seqs/combined_tree/{merge}.fa",
        ids = results + "seqs/filtering/ids/m01/{merge}.txt"
    output:
        results + "seqs/combined_tree/{merge}_m01.fa"
    shell:
        "python fasta-remove.py {input.fasta} {input.ids} {output}"

rule msa:
    input:
        results + "seqs/combined_tree/{merge}_m01.fa"
    output:
        results + "alignments/combined_tree/{merge}_m01.aln"
    shell:
        "mafft --auto --anysymbol {input} > {output}"

rule trim:
    input:
        results + "alignments/combined_tree/{merge}_m01.aln"
    output:
        results + "alignments/combined_tree/{merge}_m01.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule tree:
    """
    generate the phylogenies using iqtree
    """
    input:
        results + "alignments/combined_tree/{merge}_m01.trimal"
    output:
        results + "iqtree/combined_tree/{merge}_m01.treefile"
    params:
        uf=config["uf"],
        filename= results + "iqtree/combined_tree/{merge}_m01"
    shell:
        #"iqtree -s {input} -m MFP -mset LG -bb {params.uf} --prefix {params.filename} -redo"
        "iqtree -s {input} -m LG+R10 -bb {params.uf} --prefix {params.filename} -redo"

rule with_domain:
    input:
        results + "domain/combined_tree/{merge}.tsv",
        results + "table/phylum_table/combined_tree/{merge}_phylum.table.tsv",
        results + "iqtree/combined_tree/{merge}_m01.treefile",
        results + "seqs/combined_tree/{merge}_m01.fa"
    output:
        results + "tree_images/combined_tree/{merge}_m01.domain.tree.pdf"
    shell:
        "python ../07-tree/tree-with-domains.py {input} {output}"
