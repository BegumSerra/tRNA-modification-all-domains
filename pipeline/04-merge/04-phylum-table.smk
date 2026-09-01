# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "04-config.yaml"
results = "../../output/"
data = "../../input/"
envs = "../../envs/"
database = "../../input/database/"
PROTEIN = config["PROTEIN"]  # List of filter strings

wildcard_constraints:
    protein="|".join(PROTEIN)

rule all:
    input:
        expand(results + "table/phylum_table/{protein}_phylum.table.tsv", protein=PROTEIN)

rule phylum_table:
    input:
        fasta = results + "seqs/blast/merged/{protein}.fa",
        id_to_taxid_file = results + "table/taxa_dataset/merged_taxon_map_bact_arc.tsv",
        taxid_to_species_file = results + "table/taxa_dataset/merged_meta_data_bact_arc.tsv",
        tsv_file = results + 'table/taxa_dataset/dataset_taxa_order_updated.tsv',
        euk_taxa_file = database + "eukarya/euk_taxa.tsv"
    output:
        results + "table/phylum_table/{protein}_phylum.table.tsv"
    shell:
        "python phylum-table.py {input.fasta} {input.id_to_taxid_file} {input.taxid_to_species_file} {input.tsv_file} {output} {input.euk_taxa_file}"

