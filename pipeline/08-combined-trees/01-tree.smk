# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "config.yaml"

results = "../../output/"
data = "../../input/"

# Load configuration
MERGE = config["MERGE"]  # List of merged protein names (e.g., Tgt-QTRT1, TadA-ADAT2)
PROTEINS_MAP = config["PROTEINS"]  # Dictionary mapping merged names to list of proteins

wildcard_constraints:
    merge = "|".join(MERGE)

rule all:
    input:
        expand(results + "seqs/combined_tree/{merge}.fa", merge=MERGE),
        expand(results + "domain/combined_tree/{merge}.tsv", merge=MERGE),
        expand(results + "table/phylum_table/combined_tree/{merge}_phylum.table.tsv", merge=MERGE)
        #expand(results + "iqtree/{merge}.treefile", merge=MERGE)
        
rule merge:
    input:
        fasta = lambda wildcards: [results + f"seqs/homologs/{p}.fa" for p in PROTEINS_MAP[wildcards.merge]["fasta_names"]],
        domain = lambda wildcards: [results + "domain/blast/merged/" + p + ".tsv" for p in PROTEINS_MAP[wildcards.merge]["base_names"]],
        taxa = lambda wildcards: [results + "table/phylum_table/" + p + "_phylum.table.tsv" for p in PROTEINS_MAP[wildcards.merge]["base_names"]]
    output:
        fasta = results + "seqs/combined_tree/{merge}.fa",
        domain = results + "domain/combined_tree/{merge}.tsv",
        taxa = results + "table/phylum_table/combined_tree/{merge}_phylum.table.tsv"
    shell:
        """
        cat {input.fasta} > {output.fasta}
        cat {input.domain} > {output.domain}
        cat {input.taxa} > {output.taxa}
        """

rule msa:
    input:
        results + "seqs/combined_tree/{merge}.fa"
    output:
        results + "alignments/combined_tree/{merge}.aln"
    shell:
        "mafft --auto --anysymbol --thread {threads} {input} > {output}"

rule trim:
    input:
        results + "alignments/combined_tree/{merge}.aln"
    output:
        results + "alignments/combined_tree/{merge}.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule tree:
    """
    generate the phylogenies using iqtree
    """
    input:
        results + "alignments/combined_tree/{merge}.trimal"
    output:
        results + "iqtree/{merge}.treefile"
    params:
        uf=config["uf"],
        filename= results + "iqtree/{merge}"
    shell:
        "iqtree -s {input} -m MFP -mset LG -bb {params.uf} --prefix {params.filename} -redo"

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
