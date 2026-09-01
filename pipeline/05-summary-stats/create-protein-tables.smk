"""
Snakemake workflow to create protein-specific tables combining taxonomy, 
BLAST results, InterProScan domains, and KEGG annotations.
Snakemake workflow to highlight reference KEGG matches in protein tables

Usage:
    snakemake -s create-protein-tables.smk --cores 4
    snakemake -s create-protein-tables.smk --cores 4 --config protein=ArcS  # For specific protein
"""

import pandas as pd
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Base directories
BASE_DIR = Path("../..")
INPUT_DIR = BASE_DIR / "input"
OUTPUT_DIR = BASE_DIR / "output"
PROTEIN_TABLE_DIR = OUTPUT_DIR / "homology-analyses-tables" / "raw_tables"
logs = "logs/"

# Get list of proteins from config.yaml
PROTEINS = []
if "protein" in config:
    # Single protein specified via command line
    PROTEINS = [config["protein"]]
elif "PROTEIN" in config:
    # Use simple PROTEIN list from config file
    PROTEINS = config["PROTEIN"]
elif "PROTEIN_ARC" in config or "PROTEIN_BACT" in config or "PROTEIN_EUK" in config:
    # Use separate domain-specific protein lists from config file
    PROTEINS.extend(config.get("PROTEIN_ARC", []))
    PROTEINS.extend(config.get("PROTEIN_BACT", []))
    PROTEINS.extend(config.get("PROTEIN_EUK", []))
else:
    # Fallback: read from protein_kegg.tsv
    protein_kegg = pd.read_csv(INPUT_DIR / "kegg" / "protein_kegg.tsv")
    PROTEINS = protein_kegg['protein'].tolist()

print(f"Processing {len(PROTEINS)} proteins: {', '.join(PROTEINS)}")

# Domain types
DOMAINS = ["arc", "bact", "euk"]

# Rule to run all
rule all:
    input:
        expand(str(PROTEIN_TABLE_DIR / "{protein}_table.tsv"), protein=PROTEINS),
        expand(str(OUTPUT_DIR / "homology-analyses-tables" /"highlighted_tables" / "{protein}_highlighted.xlsx"), protein=PROTEINS)


# Rule to create protein table for a specific protein
rule create_protein_table:
    input:
        arc_metadata = INPUT_DIR / "database/archaea/archaea_gtdb_filtered_metadata.tsv",
        bact_metadata = INPUT_DIR / "database/bacteria/bacteria_gtdb_filtered_metadata.tsv",
        euk_taxa = INPUT_DIR / "database/eukarya/euk_taxa.tsv",
        arc_taxon_map = INPUT_DIR / "database/archaea/taxon_map_unique.tsv",
        bact_taxon_map = INPUT_DIR / "database/bacteria/taxon_map_unique.tsv",
        protein_list = INPUT_DIR / "kegg/protein_kegg.tsv",
        # BLAST results
        blast_arc = OUTPUT_DIR / "blast/blast/arc-arc/{protein}_unique_filtered.tab",
        blast_bact = OUTPUT_DIR / "blast/blast/bact-bact/{protein}_unique_filtered.tab",
        blast_euk = OUTPUT_DIR / "blast/blast/euk-euk/{protein}_unique_filtered.tab",
        blast_merged = OUTPUT_DIR / "blast/blast/merged/{protein}.tab",
        # Domain results
        domain_arc = OUTPUT_DIR / "domain/blast/arc-arc/{protein}_filtered_length.tsv",
        domain_bact = OUTPUT_DIR / "domain/blast/bact-bact/{protein}_filtered_length.tsv",
        domain_euk = OUTPUT_DIR / "domain/blast/euk-euk/{protein}_filtered_length.tsv",
        domain_merged = OUTPUT_DIR / "domain/blast/merged/{protein}.tsv",
        # KEGG annotations
        kegg = INPUT_DIR / "kegg/{protein}_kegg.txt"
    output:
        table = PROTEIN_TABLE_DIR / "{protein}_table.tsv"
    log:
        logs + "{protein}_table.log"
    conda:
        "pixi.toml"
    shell:
        """
        python create_protein_tables.py \
            --protein {wildcards.protein} \
            --arc-metadata {input.arc_metadata} \
            --bact-metadata {input.bact_metadata} \
            --euk-taxa {input.euk_taxa} \
            --arc-taxon-map {input.arc_taxon_map} \
            --bact-taxon-map {input.bact_taxon_map} \
            --arc-accessions {INPUT_DIR}/database/archaea/archaea.txt \
            --bact-accessions {INPUT_DIR}/database/bacteria/bacteria.txt \
            --blast-base {OUTPUT_DIR}/blast/blast \
            --domain-base {OUTPUT_DIR}/domain/blast \
            --kegg-dir {INPUT_DIR}/kegg \
            --protein-list {input.protein_list} \
            --output-dir {PROTEIN_TABLE_DIR} \
            --fasta-dir {OUTPUT_DIR}/seqs/homologs \
            > {log} 2>&1
        """

# Rule to create all protein tables at once (more efficient)
rule create_all_protein_tables:
    input:
        arc_metadata = INPUT_DIR / "database/archaea/archaea_gtdb_filtered_metadata.tsv",
        bact_metadata = INPUT_DIR / "database/bacteria/bacteria_gtdb_filtered_metadata.tsv",
        euk_taxa = INPUT_DIR / "database/eukarya/euk_taxa.tsv",
        arc_taxon_map = INPUT_DIR / "database/archaea/taxon_map_unique.tsv",
        bact_taxon_map = INPUT_DIR / "database/bacteria/taxon_map_unique.tsv",
        protein_list = INPUT_DIR / "kegg/protein_kegg.tsv",
    output:
        tables = expand(str(PROTEIN_TABLE_DIR / "{protein}_table.tsv"), protein=PROTEINS)
    log:
        logs + "all_tables.log"
    conda:
        "pixi.toml"
    params:
        proteins = " ".join(PROTEINS)
    shell:
        """
        for protein in {params.proteins}; do
            python create_protein_tables.py \
                --protein $protein \
                --arc-metadata {input.arc_metadata} \
                --bact-metadata {input.bact_metadata} \
                --euk-taxa {input.euk_taxa} \
                --arc-taxon-map {input.arc_taxon_map} \
                --bact-taxon-map {input.bact_taxon_map} \
                --arc-accessions {INPUT_DIR}/database/archaea/archaea.txt \
                --bact-accessions {INPUT_DIR}/database/bacteria/bacteria.txt \
                --blast-base {OUTPUT_DIR}/blast/blast \
                --domain-base {OUTPUT_DIR}/domain/blast \
                --kegg-dir {INPUT_DIR}/kegg \
                --protein-list {input.protein_list} \
                --output-dir {PROTEIN_TABLE_DIR} \
                --fasta-dir {OUTPUT_DIR}/seqs/homologs \
                >> {log} 2>&1
        done
        """

# Rule to generate summary statistics for all protein tables
rule summary_stats:
    input:
        tables = expand(str(PROTEIN_TABLE_DIR / "{protein}_table.tsv"), protein=PROTEINS)
    output:
        summary = PROTEIN_TABLE_DIR / "summary_statistics.tsv"
    run:
        import pandas as pd
        
        stats = []
        for table_file in input.tables:
            protein = Path(table_file).stem.replace("_table", "")
            df = pd.read_csv(table_file, sep='\t')
            
            stats.append({
                'protein': protein,
                'total_entries': len(df),
                'with_domains': df['domain'].notna().sum(),
                'with_kegg': df['kegg'].notna().sum(),
                'archaea': (df['domain_tax'] == 'Archaea').sum(),
                'bacteria': (df['domain_tax'] == 'Bacteria').sum(),
                'eukarya': (df['domain_tax'] == 'Eukarya').sum()
            })
        
        summary_df = pd.DataFrame(stats)
        summary_df.to_csv(output.summary, sep='\t', index=False)

# Rule to highlight reference KEGG matches in protein tables
rule highlight_kegg_tables:
    input:
        table = PROTEIN_TABLE_DIR / "{protein}_table.tsv",
        protein_kegg = INPUT_DIR / "kegg" / "protein_kegg.tsv"
    output:
        excel = OUTPUT_DIR / "homology-analyses-tables" / "highlighted_tables" / "{protein}_highlighted.xlsx"
    params:
        ref_dirs = [str(INPUT_DIR / "ref_seqs" / d) for d in ["arc", "bact", "euk"]]
    log:
        logs + "highlight_{protein}.log"
    conda:
        "pixi.toml"
    shell:
        """
        mkdir -p $(dirname {output.excel})
        python highlight_reference_kegg.py \
            --table {input.table} \
            --output {output.excel} \
            --ref-seqs {params.ref_dirs[0]} {params.ref_dirs[1]} {params.ref_dirs[2]} \
            --protein-kegg {input.protein_kegg} \
            > {log} 2>&1
        """