# All rights reserved

configfile: "02-config.yaml"

# Paths
results = "../../output/"
data = "../../input/"
database = "../../database/"

# Load the protein lists from the config file
F01 = config["F01"]  # List of proteins that went through f01
F00 = config.get("F00", [])  # List of proteins that did NOT go through f01

wildcard_constraints:
    protein = "|".join(F01 + F00)

rule all:
    input:
        expand(results + "seqs/filtering/fasta/{protein}_k01_f00.fa", protein=F00),
        expand(results + "seqs/filtering/fasta/{protein}_k01_f01.fa", protein=F01),
        #expand(results + "alignments/raw_alignments/{protein}_k01_f00.aln", protein=F01),
        #expand(results + "alignments/trimmed_alignments/{protein}_k01_f00.trimal", protein=F01),
        #expand(results + "alignments/raw_alignments/{protein}_k01_f01.aln", protein=F01),
        #expand(results + "alignments/trimmed_alignments/{protein}_k01_f01.trimal", protein=F01),
        "03-config.yaml"

# opy k01 filter IDs for f00 (no domain filtering for KEGG-filtered proteins)
# Only F00 proteins stop here and don't continue to domain filtering
rule copy_kegg_to_f00:
    input:
        results + "seqs/filtering/fasta/{protein}_k01.fa"
    output:
        results + "seqs/filtering/fasta/{protein}_k01_f00.fa"
    shell:
        "cp {input} {output}"

# Apply domain filtering only to F01 proteins (that already passed KEGG filtering)
# F01 proteins use the f00.fa (KEGG-filtered) FASTA and apply additional domain filtering
rule remove:
    input:
        fasta = results + "seqs/filtering/fasta/{protein}_k01.fa",
        ids = results + "seqs/filtering/ids/f01/{protein}_merged.txt"
    output:
        results + "seqs/filtering/fasta/{protein}_k01_f01.fa"
    shell:
        "python fasta-remove.py {input.fasta} {input.ids} {output}"

rule msa:
    input:
        results + "seqs/filtering/fasta/{protein}_k01_f00.fa"
    output:
        results + "alignments/raw_alignments/{protein}_k01_f00.aln"
    shell:
        "mafft --auto --thread {threads} --anysymbol {input} > {output}"

rule trim:
    input:
        results + "alignments/raw_alignments/{protein}_k01_f00.aln"
    output:
        results + "alignments/trimmed_alignments/{protein}_k01_f00.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule tree:
    """
    generate the phylogenies using iqtree
    """
    input:
        results + "alignments/trimmed_alignments/{protein}_k01_f00.trimal"
    output:
        results + "iqtree/{protein}_k01_f00.treefile"
    params:
        filename= results + "iqtree/{protein}_k01_f00"
    shell:
        #"iqtree -s {input} -m MFP -mset LG -bb {params.uf} -ntmax {threads} --prefix {params.filename} -redo"
        "iqtree -s {input} -m MFP -mset LG -bb 1000 -nt {threads} --prefix {params.filename} -redo"


rule with_domain:
    input:
        results + "domain/blast/merged/{protein}.tsv",
        results + "table/phylum_table/{protein}_phylum.table.tsv",
        results + "iqtree/{protein}_k01_f00.treefile",
        results + "seqs/filtering/fasta/{protein}_k01_f00.fa"
    output:
        results + "tree_images/{protein}_k01_f00.domain.tree.pdf"
    shell:
        "python tree-with-domains.py {input} {output}"

rule msa2:
    input:
        results + "seqs/filtering/fasta/{protein}_k01_f01.fa"
    output:
        results + "alignments/raw_alignments/{protein}_k01_f01.aln"
    shell:
        "mafft --auto --thread {threads} --anysymbol {input} > {output}"

rule trim2:
    input:
        results + "alignments/raw_alignments/{protein}_k01_f01.aln"
    output:
        results + "alignments/trimmed_alignments/{protein}_k01_f01.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule tree2:
    """
    generate the phylogenies using iqtree
    """
    input:
        results + "alignments/trimmed_alignments/{protein}_k01_f01.trimal"
    output:
        results + "iqtree/{protein}_k01_f01.treefile"
    params:
        filename= results + "iqtree/{protein}_k01_f01"
    shell:
        #"iqtree -s {input} -m MFP -mset LG -bb {params.uf} -ntmax {threads} --prefix {params.filename} -redo"
        "iqtree -s {input} -m MFP -mset LG -bb 1000 -nt {threads} --prefix {params.filename} -redo"


rule with_domain2:
    input:
        results + "domain/blast/merged/{protein}.tsv",
        results + "table/phylum_table/{protein}_phylum.table.tsv",
        results + "iqtree/{protein}_k01_f01.treefile",
        results + "seqs/filtering/fasta/{protein}_k01_f01.fa"
    output:
        results + "tree_images/{protein}_k01_f01.domain.tree.pdf"
    shell:
        "python tree-with-domains.py {input} {output}"

rule create_next_config:
    input:
        f01_files=expand(results + "seqs/filtering/fasta/{protein}_k01_f01.fa", protein=F01),
        f00_files=expand(results + "seqs/filtering/fasta/{protein}_k01_f00.fa", protein=F00)
    output:
        "03-config.yaml"
    run:
        import os
        
        # Check which proteins have e01 IDs
        e01_dir = results + "seqs/filtering/ids/e01/"
        
        f01_proteins_with_e01 = []
        f01_proteins_without_e01 = []
        for protein in F01:
            e01_file = os.path.join(e01_dir, f"{protein}.txt")
            if os.path.exists(e01_file):
                f01_proteins_with_e01.append(protein)
            else:
                f01_proteins_without_e01.append(protein)
        
        f00_proteins_with_e01 = []
        f00_proteins_without_e01 = []
        for protein in F00:
            e01_file = os.path.join(e01_dir, f"{protein}.txt")
            if os.path.exists(e01_file):
                f00_proteins_with_e01.append(protein)
            else:
                f00_proteins_without_e01.append(protein)
        
        # Write the config file
        with open(output[0], 'w') as f:
            f.write("# Configuration file for filtering eukaryotic copies from tree (step 03)\n")
            f.write("# Auto-generated based on results from step 02 (f01 filtering) and available e01 IDs\n")
            f.write("# This config defines which proteins fall into each category based on:\n")
            f.write("# 1. Proteins that completed f01 filtering (step 02)\n")
            f.write("# 2. Proteins that have e01 filtering IDs available\n\n")
            
            f.write("# Proteins that went through f01 filtering AND have e01 IDs\n")
            f.write("# These will be processed as {protein}_f01 -> {protein}_f01_e01\n")
            f.write("# Source: Proteins with both fasta/*_f01.fa files and ids/e01/{protein}.txt files\n")
            f.write("K01_F01_E01:\n")
            if f01_proteins_with_e01:
                for p in f01_proteins_with_e01:
                    f.write(f"  - {p}\n")
            else:
                f.write("  []\n")
            
            f.write("\n# Proteins that did NOT go through f01 filtering (f00) AND have e01 IDs\n")
            f.write("# These will be processed as {protein}_f00 -> {protein}_f00_e01\n")
            f.write("# Source: Proteins with both fasta/*_f00.fa files and ids/e01/{protein}.txt files\n")
            f.write("K01_F00_E01:\n")
            if f00_proteins_with_e01:
                for p in f00_proteins_with_e01:
                    f.write(f"  - {p}\n")
            else:
                f.write("  []\n")
            
            
            f.write("\n# Proteins that went through f01 filtering but do NOT have e01 IDs\n")
            f.write("# These will be copied as {protein}_f01 -> {protein}_f01_e00\n")
            f.write("# Source: Proteins with fasta/*_f01.fa files but no ids/e01/{protein}.txt files\n")
            f.write("K01_F01_E00:\n")
            if f01_proteins_without_e01:
                for p in f01_proteins_without_e01:
                    f.write(f"  - {p}\n")
            else:
                f.write("  []\n")
            
            f.write("\n# Proteins that did NOT go through f01 filtering (f00) and do NOT have e01 IDs\n")
            f.write("# These will be copied as {protein}_f00 -> {protein}_f00_e00\n")
            f.write("# Source: Proteins with fasta/*_f00.fa files but no ids/e01/{protein}.txt files\n")
            f.write("K01_F00_E00:\n")
            if f00_proteins_without_e01:
                for p in f00_proteins_without_e01:
                    f.write(f"  - {p}\n")
            else:
                f.write("  []\n")
            
            f.write("\n# IQTree ultrafast bootstrap parameter\n")
            f.write("uf: 1000\n")

