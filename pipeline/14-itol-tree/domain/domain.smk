# Snakefile for creating iTOL domain tables from InterProScan output
# Copyright © 2024, Begüm Serra Büyüktarakçı

import os

configfile: "config.yaml"

PROTEIN = config["PROTEIN"]
COMBINED_PROTEIN = config.get("COMBINED_PROTEIN", [])
protein_to_fasta = config["MAP"]
protein_to_domain = config.get("DOMAIN_MAP", {})

# If DOMAIN_MAP not provided, use same name as protein
for p in PROTEIN + COMBINED_PROTEIN:
    if p not in protein_to_domain:
        protein_to_domain[p] = p

wildcard_constraints:
    protein = "|".join(PROTEIN),
    combined = "|".join(COMBINED_PROTEIN) if COMBINED_PROTEIN else "^$"

def get_fasta(wildcards):
    """Get FASTA file for a protein"""
    return f"{config['fasta']}/{protein_to_fasta[wildcards.protein]}.fa"

def get_fasta_combined(wildcards):
    """Get FASTA file for a combined protein"""
    return f"{config['fasta_combined']}/{protein_to_fasta[wildcards.combined]}.fa"

def get_domain_tsv(wildcards):
    """Get InterProScan TSV file for a protein"""
    domain_name = protein_to_domain.get(wildcards.protein, wildcards.protein)
    return f"{config['domain']}/{domain_name}.tsv"

def get_domain_tsv_combined(wildcards):
    """Get InterProScan TSV file for a combined protein"""
    domain_name = protein_to_domain.get(wildcards.combined, wildcards.combined)
    return f"{config['domain']}/{domain_name}.tsv"

rule all:
    input:
        expand(config["output"] + "/{protein}_domains.tsv", protein=PROTEIN),
        expand(config["output"] + "/combined_tree/{combined}_domains.tsv", combined=COMBINED_PROTEIN),
        expand(config["output"] + "/{protein}_domains_itol.txt", protein=PROTEIN),
        expand(config["output"] + "/combined_tree/{combined}_domains_itol.txt", combined=COMBINED_PROTEIN)

rule create_domain_table:
    input:
        fasta = get_fasta,
        domain = get_domain_tsv
    output:
        table = config["output"] + "/{protein}_domains.tsv"
    params:
        protein = lambda wildcards: wildcards.protein
    shell:
        """
        python create_domain_table.py {input.domain} \
            -f {input.fasta} \
            -o {output.table} \
            --legend
        """

rule create_domain_table_combined:
    input:
        fasta = get_fasta_combined,
        domain = get_domain_tsv_combined
    output:
        table = config["output"] + "/combined_tree/{combined}_domains.tsv"
    params:
        combined = lambda wildcards: wildcards.combined
    shell:
        """
        python create_domain_table.py {input.domain} \
            -f {input.fasta} \
            -o {output.table} \
            --legend
        """

rule convert_to_itol:
    input:
        tsv = config["output"] + "/{protein}_domains.tsv"
    output:
        itol = config["output"] + "/{protein}_domains_itol.txt"
    params:
        label = lambda wildcards: wildcards.protein
    shell:
        """
        python tsv_to_itol_domains.py {input.tsv} \
            -o {output.itol} \
            -l {params.label}
        """

rule convert_to_itol_combined:
    input:
        tsv = config["output"] + "/combined_tree/{combined}_domains.tsv"
    output:
        itol = config["output"] + "/combined_tree/{combined}_domains_itol.txt"
    params:
        label = lambda wildcards: wildcards.combined
    shell:
        """
        python tsv_to_itol_domains.py {input.tsv} \
            -o {output.itol} \
            -l {params.label}
        """