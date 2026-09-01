# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

"""
Add archaeal sequences based on a01 ID list.
This adds sequences with kegg=False but e-value below specified threshold.
- Tuc1: e-value < 1e-25
- TRM9B: e-value < 1e-20
"""

configfile: "add.yaml"

results = "../../output/"
database = "../../input/database/"

# Define proteins and their e-value thresholds
PROTEIN_EVALUE = {
    "Tuc1": 1e-25,
    "TRM9B": 1e-20  # Less strict threshold for TRM9B
}

# Define which proteins should include all domains (not just archaea)
ALL_DOMAINS = ["TRM9B"]

rule all:
    input:
        expand(results + "seqs/homologs/{protein}_k01_f01_e00_c01_nogaps_r01_a01.fa", protein=PROTEIN_EVALUE.keys())

rule extract_ids:
    """
    Extract IDs from highlighted Excel file.
    Tuc1: Filter for archaeal sequences with kegg=False and e-value < 1e-25.
    TRM9B: Filter for all domain sequences with kegg=False and e-value < 1e-20.
    """
    input:
        excel = results + "homology-analyses-tables/highlighted_tables/{protein}_highlighted.xlsx"
    output:
        id_list = results + "seqs/filtering/ids/a01/{protein}.txt"
    params:
        evalue = lambda wildcards: PROTEIN_EVALUE[wildcards.protein],
        domain_flag = lambda wildcards: "all" if wildcards.protein in ALL_DOMAINS else "archaea"
    shell:
        """
        python find-ids.py {input.excel} {output.id_list} {params.evalue} {params.domain_flag}
        """

rule add_sequences:
    """
    Add archaeal sequences from a01 ID list to the fasta file.
    These are sequences with kegg=False but good e-values (< 1e-25).
    """
    input:
        fasta = results + "seqs/homologs/Tuc1_k01_f01_e00_c01_nogaps_r01.fa",
        id_list = results + "seqs/filtering/ids/a01/Tuc1.txt"
    output:
        fasta_out = results + "seqs/homologs/Tuc1_k01_f01_e00_c01_nogaps_r01_a01.fa"
    params:
        db_dir = database
    shell:
        """
        python add-seqs.py {input.fasta} {input.id_list} {output.fasta_out} {params.db_dir}
        """

rule extract_sequences:
    """
    Add sequences for TRM9B from a01 ID list to the fasta file.
    This includes all domain sequences with kegg=False but good e-values (< 1e-20).
    """
    input:
        fasta = results + "seqs/homologs/TRM9B_k01_f01_e00_c01_nogaps_r01.fa",
        id_list = results + "seqs/filtering/ids/a01/TRM9B.txt"
    output:
        fasta_out = results + "seqs/homologs/TRM9B_k01_f01_e00_c01_nogaps_r01_a01.fa"
    params:
        db_dir = database
    shell:
        """
        python add-seqs.py {input.fasta} {input.id_list} {output.fasta_out} {params.db_dir}
        """