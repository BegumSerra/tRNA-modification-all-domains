# Copyright (c) 2024, Begüm Serra Büyüktarakçı
# All rights reserved.

results = "../../output/"
data = "../../input/"
envs = "../../envs/"
database = "../../input/database/"

rule all:
    input: 
        results + "table/taxa_dataset/merged_meta_data_bact_arc.tsv",
        results + "table/taxa_dataset/merged_taxon_map_bact_arc.tsv",
        results + "table/taxa_dataset/dataset_taxa_order_updated.tsv",
        results + "table/taxa_dataset/color_preview.html",
        results + "table/protein_genome_mapping.tsv"

rule merge_meta_data_bacteria_archaea:
    input:
        database + "bacteria/bacteria_gtdb_filtered_metadata.tsv",
        database + "archaea/archaea_gtdb_filtered_metadata.tsv"
    output:
        results + "table/taxa_dataset/merged_meta_data_bact_arc.tsv"
    shell:
        """
        python merge_metadata.py --bacteria {input[0]} --archaea {input[1]} --output {output}
        """

rule merge_taxon_map_bacteria_archaea:
    input:
        database + "bacteria/taxon_map_unique.tsv",
        database + "archaea/taxon_map_unique.tsv"
    output:
        results + "table/taxa_dataset/merged_taxon_map_bact_arc.tsv"
    shell:
        """
        cat {input[0]} {input[1]} > {output}
        """

rule update_taxa_order:
    input:
        results + "table/taxa_dataset/merged_meta_data_bact_arc.tsv",
        data + "table/dataset_taxa_order_manual.txt"
    output:
        results + "table/taxa_dataset/dataset_taxa_order_updated.tsv"
    shell:
        """
        python update_taxonomy.py {input} {output}
        """

rule color_taxa_preview:
    input:
        results + "table/taxa_dataset/dataset_taxa_order_updated.tsv",
        data + "table/color.txt"
    output:
        results + "table/taxa_dataset/color_preview.html"
    shell:
        """
        python create-color-file.py {input[0]} {input[1]} {output}
        """

rule mapping_protein_genome:
    input:
        results + "table/taxa_dataset/merged_meta_data_bact_arc.tsv"
    output:
        results + "table/protein_genome_mapping.tsv"
    shell:
        """
        python create_protein_genome_mapping.py {database} {input} {output}
        """