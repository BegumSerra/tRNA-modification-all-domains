# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

results = "../../../output/"
data = "../../../input/"

rule all:
    input:
        results + "table/presence_absence/dot_graph/main/circle-presence-function.pdf"

rule reorder:
    input:
        results + "table/presence_absence/protein_presence_ordered_by_function.tsv"
    output:
        results + "table/presence_absence/dot_graph/main/protein_presence_ordered_by_group.xlsx"
    shell:
        "python ../suppl-dot/reorder-columns.py {input} {output}"

rule binary_data:
    input:
        xlsx = results + "table/presence_absence/dot_graph/main/protein_presence_ordered_by_group.xlsx",
        order = "order.txt"
    output:
        results + "table/presence_absence/dot_graph/main/protein_presence_binary.csv"
    shell:
        "python ../suppl-dot/binary-data.py {input.xlsx} {input.order} {output}"

rule dot_graph:
    input:
        csv = results + "table/presence_absence/dot_graph/main/protein_presence_binary.csv",
        order = "order.txt"
    output:
        results + "table/presence_absence/dot_graph/main/circle-presence-function.pdf"
    shell:
        "python circle-presence-function.py {input.csv} {input.order} {output}"

