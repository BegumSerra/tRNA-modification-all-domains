# Copyright (c) 2024, Begüm Serra Büyüktarakçı
# All rights reserved.

"""
Snakemake workflow to download protein sequences and metadata from NCBI
for bacteria and archaea assemblies
"""

from pathlib import Path

# Configuration
DOMAINS = ["bacteria", "archaea"]

# Get the project root (two levels up from this file)
WORKFLOW_DIR = Path(workflow.basedir)
PROJECT_ROOT = WORKFLOW_DIR.parent.parent

# Target rule - runs everything
rule all:
    input:
        # Combined protein files
        expand(str(PROJECT_ROOT / "input/database/{domain}/{domain}_ncbi_proteins.faa"), domain=DOMAINS),
        # Deduplicated protein files
        expand(str(PROJECT_ROOT / "input/database/{domain}/{domain}_ncbi_proteins_unique.faa"), domain=DOMAINS),
        # GTDB metadata files
        expand(str(PROJECT_ROOT / "input/database/{domain}/gtdb_{domain}_metadata.tsv"), domain=DOMAINS),
        # Filtered GTDB metadata
        expand(str(PROJECT_ROOT / "input/database/{domain}/{domain}_gtdb_filtered_metadata.tsv"), domain=DOMAINS),
        # Taxon map files
        expand(str(PROJECT_ROOT / "input/database/{domain}/taxon_map_unique.tsv"), domain=DOMAINS)

# Download protein sequences from NCBI
rule download_proteins:
    input:
        accessions=str(PROJECT_ROOT / "input/database/{domain}/{domain}.txt")
    output:
        proteins=str(PROJECT_ROOT / "input/database/{domain}/{domain}_ncbi_proteins.faa")
    params:
        download_dir=str(PROJECT_ROOT / "input/database/{domain}/{domain}_downloads"),
        batch_size=50
    log:
        str("logs/download_proteins_{domain}.log")
    shell:
        """
        python download_ncbi_proteins.py \
            {input.accessions} \
            {params.download_dir} \
            {output.proteins} \
            --batch-size {params.batch_size} \
            2>&1 | tee {log}
        """

# Download GTDB taxonomy metadata
rule download_gtdb_taxonomy:
    output:
        gtdb_metadata=str(PROJECT_ROOT / "input/database/{domain}/gtdb_{domain}_metadata.tsv")
    params:
        domain="{domain}"
    log:
        str("logs/download_gtdb_{domain}.log")
    shell:
        """
        python download_gtdb_taxonomy.py \
            --domain {params.domain} \
            --output {output.gtdb_metadata} \
            2>&1 | tee {log}
        """

# Filter GTDB metadata to only our assemblies
rule filter_gtdb_metadata:
    input:
        gtdb_metadata=str(PROJECT_ROOT / "input/database/{domain}/gtdb_{domain}_metadata.tsv"),
        accession_list=str(PROJECT_ROOT / "input/database/{domain}/{domain}.txt")
    output:
        filtered_metadata=str(PROJECT_ROOT / "input/database/{domain}/{domain}_gtdb_filtered_metadata.tsv")
    log:
        str("logs/filter_gtdb_{domain}.log")
    shell:
        """
        python filter_gtdb_metadata.py \
            --gtdb-metadata {input.gtdb_metadata} \
            --accession-list {input.accession_list} \
            --output {output.filtered_metadata} \
            2>&1 | tee {log}
        
        # Clean RS_/GB_ prefixes from accession column
        python -c "
import pandas as pd
df = pd.read_csv('{output.filtered_metadata}', sep='\\t')
df['accession'] = df['accession'].str.replace('RS_', '').str.replace('GB_', '')
df.to_csv('{output.filtered_metadata}', sep='\\t', index=False)
print(f'Cleaned {{len(df)}} assemblies')
        " 2>&1 | tee -a {log}
        """

# Create taxon_map.tsv mapping protein IDs to taxonomy IDs
rule create_taxon_map:
    input:
        filtered_metadata=str(PROJECT_ROOT / "input/database/{domain}/{domain}_gtdb_filtered_metadata.tsv")
    output:
        taxon_map=str(PROJECT_ROOT / "input/database/{domain}/taxon_map.tsv")
    params:
        download_dir=lambda wildcards: str(PROJECT_ROOT / f"input/database/{wildcards.domain}/ncbi_dataset_{'130824' if wildcards.domain == 'bacteria' else '190824'}")
    log:
        str("logs/create_taxon_map_{domain}.log")
    shell:
        """
        python create_taxon_map.py \
            {params.download_dir} \
            {input.filtered_metadata} \
            {output.taxon_map} \
            2>&1 | tee {log}
        """

# Deduplicate proteins and update taxon_map
rule deduplicate_proteins:
    input:
        proteins=str(PROJECT_ROOT / "input/database/{domain}/{domain}_ncbi_proteins.faa"),
        taxon_map=str(PROJECT_ROOT / "input/database/{domain}/taxon_map.tsv")
    output:
        proteins_unique=str(PROJECT_ROOT / "input/database/{domain}/{domain}_ncbi_proteins_unique.faa"),
        taxon_map_unique=str(PROJECT_ROOT / "input/database/{domain}/taxon_map_unique.tsv")
    log:
        str("logs/deduplicate_{domain}.log")
    shell:
        """
        python deduplicate_proteins.py \
            --input-fasta {input.proteins} \
            --input-taxon-map {input.taxon_map} \
            --output-fasta {output.proteins_unique} \
            --output-taxon-map {output.taxon_map_unique} \
            2>&1 | tee {log}
        """
