# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

results = "../../../output/"
data = "../../../input/"

rule all:
    input:
        results + "table/presence_absence/dot_graph/main/circle_presence_function_new_grouping_filtered_proteins2.pdf"

rule binary_data:
    input:
        xlsx = results + "table/presence_absence/new_grouping/colored_presence_only.xlsx",
    output:
        results + "table/presence_absence/dot_graph/main/protein_presence_binary_new_grouping.csv"
    shell:
        "python binary-data-new-grouping.py {input.xlsx} {output}"

rule dot_graph:
    input:
        csv = results + "table/presence_absence/dot_graph/main/protein_presence_binary_new_grouping.csv",
        groups = "group_function2.txt"
    output:
        results + "table/presence_absence/dot_graph/main/circle_presence_function_new_grouping_filtered_proteins2.pdf"
    shell:
        "python circle-presence-function-new-grouping.py {input.csv} {input.groups} {output}"

