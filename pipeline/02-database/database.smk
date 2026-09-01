# Copyright (c) 2024, Begüm Serra Büyüktarakçı
# All rights reserved.

results = "../../output/"
data = "../../input/"
envs = "../../envs/"
database = "../../input/database/"
logs = "logs/"

rule all:
    input:
        results + "blastdb/archaea/blastdb.done",
        results + "blastdb/bacteria/blastdb.done",
        results + "blastdb/eukarya/blastdb.done"


rule makeblastdb_archaea:
    input:
        fasta = database + "archaea/archaea_ncbi_proteins_unique.faa",
        taxon = database + "archaea/taxon_map_unique.tsv"
    output:
        done = results + "blastdb/archaea/blastdb.done"
    log:
        logs + "makeblastdb_archaea.log"
    shell:
        """
        mkdir -p {results}blastdb/archaea/
        makeblastdb -in {input.fasta} -dbtype prot -parse_seqids -taxid_map {input.taxon} -out {results}blastdb/archaea/blastdb 2>&1 | tee {log}
        touch {output.done}
        """

rule makeblastdb_bacteria:
    input:
        fasta = database + "bacteria/bacteria_ncbi_proteins_unique.faa",
        taxon = database + "bacteria/taxon_map_unique.tsv"
    output:
        done = results + "blastdb/bacteria/blastdb.done"
    log:
        logs + "makeblastdb_bacteria.log"
    shell:
        """
        mkdir -p {results}blastdb/bacteria/
        makeblastdb -in {input.fasta} -dbtype prot -parse_seqids -taxid_map {input.taxon} -out {results}blastdb/bacteria/blastdb 2>&1 | tee {log}
        touch {output.done}
        """

rule makeblastdb_eukarya:
    input:
        fasta = database + "eukarya/euk_proteomes.faa"
    output:
        done = results + "blastdb/eukarya/blastdb.done"
    log:
        logs + "makeblastdb_eukarya.log"
    shell:
        """
        mkdir -p {results}blastdb/eukarya/
        makeblastdb -in {input.fasta} -dbtype prot -out {results}blastdb/eukarya/blastdb 2>&1 | tee {log}
        touch {output.done}
        """
