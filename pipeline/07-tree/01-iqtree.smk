# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "tree-TgtA.yaml"
wildcard_constraints:
    protein="|".join(config["PROTEIN"])
results = "../../output/"
data = "../../input/"
envs = "../../envs/"
database = "../../input/database/"
PROTEIN = config["PROTEIN"]
filter_to_protein = config["FILTER_TO_PROTEIN"]

rule all:
    input:
        expand(results + "iqtree/{protein}_nogaps.treefile", protein=PROTEIN)

rule msa:
    input:
        results + "seqs/filtering/fasta/{protein}.fa"
    output:
        results + "alignments/raw_alignments/{protein}.aln"
    threads:
        10
    shell:
        "mafft --auto --anysymbol --thread {threads} {input} > {output}"

rule trim:
    input:
        results + "alignments/raw_alignments/{protein}.aln"
    output:
        results + "alignments/trimmed_alignments/{protein}.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule filter_gaps_trimal:
    input:
        align = results + "alignments/trimmed_alignments/{protein}.trimal",
        seqs = results + "seqs/filtering/fasta/{protein}.fa"
    output:
        results + "seqs/homologs/{protein}_nogaps.fa"

    shell:
        "python filter-gaps-trimal.py {input.align} {input.seqs} {output}"

rule msa2:
    input:
        results + "seqs/homologs/{protein}_nogaps.fa"
    output:
        results + "alignments/raw_alignments/{protein}_nogaps.aln"
    threads:
        10
    shell:
        "mafft --auto --anysymbol --thread {threads} {input} > {output}"

rule trim2:
    input:
        results + "alignments/raw_alignments/{protein}_nogaps.aln"
    output:
        results + "alignments/trimmed_alignments/{protein}_nogaps.trimal"
    threads:
        1
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule tree:
    """
    generate the phylogenies using iqtree
    """
    input:
        results + "alignments/trimmed_alignments/{protein}_nogaps.trimal"
    output:
        results + "iqtree/{protein}_nogaps.treefile"
    params:
        uf=config["uf"],
        filename= results + "iqtree/{protein}_nogaps"
    shell:
        "iqtree -s {input} -m MFP -mset LG -bb {params.uf} -nt AUTO --prefix {params.filename} -redo"

