# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "config.yaml"
wildcard_constraints:
    protein = "|".join(config["PROTEIN"])

rule all:
    input:
         #expand(config["range"] +"/{protein}_dataset_ranges_deeploc_C_N_M_P.txt", protein=config["PROTEIN"]),
         #expand(config["range"] +"/{protein}_dataset_ranges_deeploc_CN_M_P.txt", protein=config["PROTEIN"]),
         #expand(config["range"] +"/{protein}_dataset_ranges_deeploc_CN_MP.txt", protein=config["PROTEIN"]),
         expand(config["range"] +"/{protein}_dataset_ranges_deeploc_threshold.txt", protein=config["PROTEIN"])

# rule deeploc_filter_C_N_M_P:
#     input:
#         fasta = config["table"] +"/{protein}.csv"
#     output:
#         filtered = config["range"] +"/{protein}_dataset_ranges_deeploc_C_N_M_P.txt"
#     shell:
#         "python filter-deeploc-C-N-M-P.py {input.fasta} {output.filtered}"

# rule deeploc_filter_CN_M_P:
#     input:
#         fasta = config["table"] +"/{protein}.csv"
#     output:
#         filtered = config["range"] +"/{protein}_dataset_ranges_deeploc_CN_M_P.txt"
#     shell:
#         "python filter-deeploc-combined-CN-M-P.py {input.fasta} {output.filtered}"

# rule deeploc_filter_CN_MP:
#     input:
#         fasta = config["table"] +"/{protein}.csv"
#     output:
#         filtered = config["range"] +"/{protein}_dataset_ranges_deeploc_CN_MP.txt"
#     shell:
#         "python filter-deeploc-combined-CN-MP.py {input.fasta} {output.filtered}"


rule deeploc_filter_threshold:
    input:
        fasta = config["table"] +"/{protein}_deeploc.csv"
    output:
        filtered = config["range"] +"/{protein}_dataset_ranges_deeploc_threshold.txt"
    shell:
        "python filter-deeploc-threshold.py {input.fasta} {output.filtered}"
