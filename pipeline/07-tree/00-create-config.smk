# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

results = "../../output/"
fasta_dir = results + "seqs/filtering/fasta"

rule all:
    input:
        "tree.yaml"

rule generate_config:
    input:
        fasta_dir
    output:
        "tree.yaml"
    shell:
        "python generate_tree_config.py {input} {output}"
