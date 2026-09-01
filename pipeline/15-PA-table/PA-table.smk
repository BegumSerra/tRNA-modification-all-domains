# Copyright © 2026, Begüm Serra Büyüktarakçı
# All rights reserved

results = "../../output/"
data = "../../input/"
envs = "../../envs/"
database = "../../database/"


rule all:
    input:
        results + "table/presence_absence/protein_presence_colored.xlsx",
        results + "table/presence_absence/new_grouping/colored_all_presence_absence.xlsx",
        results + "table/presence_absence/new_grouping/colored_presence_only.xlsx"

rule create_table:
    input:
        results + "table/taxa_dataset/dataset_taxa_order_updated.tsv",
        results + "table/protein_genome_mapping.tsv"
    params:
        results + "seqs/homologs",
        results + "seqs/combined_tree"
    output:
        results + "table/presence_absence/dataset_with_protein_presence.tsv",
        results + "table/presence_absence/protein_presence_statistics.tsv"
    log:
        results + "table/presence_absence/create_table.log"
    shell:
        "python analyse_protein_presence.py {input[0]} {params[0]} {params[1]} {input[1]} {output[0]} {output[1]} > {log} 2>&1"

rule order_proteins:
    input:
        data = results + "table/presence_absence/dataset_with_protein_presence.tsv",
        order = "protein_order_function.txt"
    output:
        results + "table/presence_absence/protein_presence_ordered_by_function.tsv"
    shell:
        "python order-proteins.py {input.data} {input.order} {output}"

rule color_table:
    input:
        data = results + "table/presence_absence/protein_presence_ordered_by_function.tsv",
        order = data + "table/color.txt"
    output:
        results + "table/presence_absence/protein_presence_colored.xlsx"
    shell:
        "python color-PA.py {input.data} {input.order} {output}"

rule color_table_by_group:
    input:
        data = results + "table/presence_absence/protein_presence_ordered_by_function.tsv",
        order = results + "table/taxa_dataset/taxa_color_manual_order.tsv"
    output:
        results + "table/presence_absence/new_grouping/colored_all_presence_absence.xlsx"
    shell:
        "python order-color-PA-by-group.py {input.data} {input.order} {output}"

rule color_table_ones_only:
    input:
        data = results + "table/presence_absence/protein_presence_ordered_by_function.tsv",
        order = results + "table/taxa_dataset/taxa_color_manual_order.tsv"
    output:
        results + "table/presence_absence/new_grouping/colored_presence_only.xlsx"
    shell:
        "python color-PA-ones-only.py {input.data} {input.order} {output}"

        
