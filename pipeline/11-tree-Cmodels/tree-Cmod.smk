# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

import pandas as pd

configfile: "tree.yaml"

results = "../../output/"
data = "../../input/"
envs = "../../envs/"
database = "../../input/database/"
PROTEIN = config["PROTEIN"]
filter_to_protein = config["FILTER_TO_PROTEIN"]

# Separate combined and regular proteins
COMBINED_PROTEINS = [p for p in PROTEIN if '-' in p]
REGULAR_PROTEINS = [p for p in PROTEIN if '-' not in p]

wildcard_constraints:
    protein="|".join(REGULAR_PROTEINS) if REGULAR_PROTEINS else "^$",
    combined="|".join(COMBINED_PROTEINS) if COMBINED_PROTEINS else "^$"

# Function to read models from tree_models_summary.tsv
def get_best_model(protein_name, is_combined=False):
    """
    Get the best model for a protein from tree_models_summary.tsv
    """
    try:
        # Read the summary file
        df = pd.read_csv("tree_models_summary.tsv", sep='\t')
        
        if is_combined:
            # For combined proteins, the protein_name already includes the full suffix
            # e.g., protein_name = "Tgt-QTRT1_m01_r01"
            row = df[df['protein'] == protein_name]
        else:
            # For regular proteins, look for pattern with _k01_f*_e*_c*_nogaps_r01 suffix
            matching_rows = df[df['protein'].str.contains(f"^{protein_name}_k01_", na=False)]
            row = matching_rows.head(1)  # Take first match if multiple
        
        if not row.empty and pd.notna(row.iloc[0]['best_model']):
            return row.iloc[0]['best_model']
        else:
            print(f"Warning: No model found for {protein_name}, using LG+R7 as default")
            return "LG+R7"
    except Exception as e:
        print(f"Error reading tree models summary for {protein_name}: {e}")
        return "LG+R7"

rule all:
    input:
        expand(results + "iqtree/{protein}_Cmod.treefile", protein=REGULAR_PROTEINS),
        expand(results + "iqtree/combined_tree/{combined}_Cmod.treefile", combined=COMBINED_PROTEINS)

rule msa:
    input:
        results + "seqs/homologs/{protein}.fa"
    output:
        results + "alignments/raw_alignments/{protein}.aln"
    threads:
        10
    shell:
        "mafft --auto --anysymbol --thread {threads} {input} > {output}"

rule msa_combined:
    input:
        results + "seqs/combined_tree/{combined}.fa"
    output:
        results + "alignments/raw_alignments/{combined}.aln"
    threads:
        10
    shell:
        "mafft --auto --anysymbol --thread {threads} {input} > {output}"

rule trim:
    input:
        results + "alignments/raw_alignments/{protein}.aln"
    output:
        results + "alignments/trimmed_alignments/{protein}.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule trim_combined:
    input:
        results + "alignments/raw_alignments/{combined}.aln"
    output:
        results + "alignments/trimmed_alignments/{combined}.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule tree:
    """
    Generate the phylogenies using iqtree for regular proteins with pre-determined best model.
    """
    input:
        results + "alignments/trimmed_alignments/{protein}.trimal"
    output:
        results + "iqtree/{protein}_Cmod.treefile"
    params:
        uf=1000,
        filename=results + "iqtree/{protein}_Cmod",
        model=lambda wildcards: get_best_model(wildcards.protein, is_combined=False)
    shell:
        "iqtree -s {input} -m {params.model} -bb {params.uf} -nt AUTO --prefix {params.filename} -redo"

rule tree_combined:
    """
    Generate the phylogenies using iqtree for combined tree proteins with pre-determined best model.
    """
    input:
        results + "alignments/trimmed_alignments/{combined}.trimal"
    output:
        results + "iqtree/combined_tree/{combined}_Cmod.treefile"
    params:
        uf=1000,
        filename=results + "iqtree/combined_tree/{combined}_Cmod",
        model=lambda wildcards: get_best_model(wildcards.combined, is_combined=True)
    shell:
        "iqtree -s {input} -m {params.model} -bb {params.uf} -nt AUTO --prefix {params.filename} -redo"

