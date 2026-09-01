# All rights reserved

configfile: "03-config.yaml"

# Paths
results = "../../output/"
data = "../../input/"
database = "../../database/"

# Constrain wildcards to prevent greedy matching
# protein should not contain underscores (protein names like ADAT2, not ADAT2_k01)
# type should match k01_f01 or k01_f00 patterns
wildcard_constraints:
    protein="[^_]+",  # No underscores in protein name
    type="k01_f(01|00)"  # Must be k01_f01 or k01_f00

# Create lists with proper type annotations
# Input files from step 02 use "k01_f01" and "k01_f00" naming
K01_F01_E01_list = [(p, "k01_f01") for p in config["K01_F01_E01"]]
K01_F00_E01_list = [(p, "k01_f00") for p in config["K01_F00_E01"]]
K01_F01_E00_list = [(p, "k01_f01") for p in config["K01_F01_E00"]]
K01_F00_E00_list = [(p, "k01_f00") for p in config["K01_F00_E00"]]

rule all:
    input:
        expand(results + "seqs/filtering/fasta/{protein}_{type}_e01.fa",
               zip, protein=[p for p, t in K01_F01_E01_list], type=[t for p, t in K01_F01_E01_list]),
        expand(results + "seqs/filtering/fasta/{protein}_{type}_e01.fa",
               zip, protein=[p for p, t in K01_F00_E01_list], type=[t for p, t in K01_F00_E01_list]),
        expand(results + "seqs/filtering/fasta/{protein}_{type}_e00.fa",
               zip, protein=[p for p, t in K01_F01_E00_list], type=[t for p, t in K01_F01_E00_list]),
        expand(results + "seqs/filtering/fasta/{protein}_{type}_e00.fa",
               zip, protein=[p for p, t in K01_F00_E00_list], type=[t for p, t in K01_F00_E00_list]),
        #expand(results + "alignments/raw_alignments/{protein}_{type}_e01.aln",
        #       zip, protein=[p for p, t in K01_F01_E01_list], type=[t for p, t in K01_F01_E01_list]),
        #expand(results + "alignments/raw_alignments/{protein}_{type}_e01.aln",
        #       zip, protein=[p for p, t in K01_F00_E01_list], type=[t for p, t in K01_F00_E01_list]),
        #expand(results + "alignments/trimmed_alignments/{protein}_{type}_e01.trimal",
        #       zip, protein=[p for p, t in K01_F01_E01_list], type=[t for p, t in K01_F01_E01_list]),
        #expand(results + "alignments/trimmed_alignments/{protein}_{type}_e01.trimal",
        #       zip, protein=[p for p, t in K01_F00_E01_list], type=[t for p, t in K01_F00_E01_list]),
        "04-config.yaml"

rule copy_no_e01:
    input:
        lambda wildcards: results + f"seqs/filtering/fasta/{wildcards.protein}_{wildcards.type}.fa"
    output:
        results + "seqs/filtering/fasta/{protein}_{type}_e00.fa"
    shell:
        "cp {input} {output}"

rule remove:
    input:
        fasta=lambda wildcards: results + f"seqs/filtering/fasta/{wildcards.protein}_{wildcards.type}.fa",
        ids=results + "seqs/filtering/ids/e01/{protein}.txt"
    output:
        results + "seqs/filtering/fasta/{protein}_{type}_e01.fa"
    shell:
        "python fasta-remove.py {input.fasta} {input.ids} {output} || {{ echo 'Error processing {wildcards.protein}'; exit 1; }}"

rule msa:
    input:
        results + "seqs/filtering/fasta/{protein}_{type}_e01.fa"
    output:
        results + "alignments/raw_alignments/{protein}_{type}_e01.aln"
    shell:
        "mafft --auto --anysymbol --thread {threads} {input} > {output}"

rule trim:
    input:
        results + "alignments/raw_alignments/{protein}_{type}_e01.aln"
    output:
        results + "alignments/trimmed_alignments/{protein}_{type}_e01.trimal"
    shell:
        "trimal -in {input} -out {output} -fasta -gt 0.5"

rule tree:
    """
    generate the phylogenies using iqtree
    """
    input:
        results + "alignments/trimmed_alignments/{protein}_{type}_e01.trimal"
    output:
        results + "iqtree/{protein}_{type}_e01.treefile"
    params:
        uf=config["uf"],
        filename= results + "iqtree/{protein}_{type}_e01"
    shell:
        #"iqtree -s {input} -m MFP -mset LG -bb {params.uf} -nt {threads} -ntmax {threads} --prefix {params.filename} -redo"
        "iqtree -s {input} -m MFP -mset LG -bb {params.uf} --prefix {params.filename} -redo"

rule with_domain:
    input:
        results + "domain/blast/merged/{protein}.tsv",
        results + "table/phylum_table/{protein}_phylum.table.tsv",
        results + "iqtree/{protein}_{type}_e01.treefile",
        results + "seqs/filtering/fasta/{protein}_{type}_e01.fa"
    output:
        results + "tree_images/{protein}_{type}_e01.domain.tree.pdf"
    shell:
        "python tree-with-domains.py {input} {output}"

rule create_next_config:
    input:
        e01_files=expand(results + "seqs/filtering/fasta/{protein}_{type}_e01.fa",
                        zip, protein=config["K01_F01_E01"] + config["K01_F00_E01"],
                        type=["k01_f01"]*len(config["K01_F01_E01"]) + ["k01_f00"]*len(config["K01_F00_E01"])),
        e00_files=expand(results + "seqs/filtering/fasta/{protein}_{type}_e00.fa",
                        zip, protein=config["K01_F01_E00"] + config["K01_F00_E00"],
                        type=["k01_f01"]*len(config["K01_F01_E00"]) + ["k01_f00"]*len(config["K01_F00_E00"]))
    output:
        "04-config.yaml"
    run:
        import os
        import yaml
        
        # Collect all proteins and their filters
        protein_list = []
        filter_list = []
        protein_to_filter = {}
        
        # Add e01 trees
        for protein, ftype in zip(config["K01_F01_E01"] + config["K01_F00_E01"],
                                   ["k01_f01"]*len(config["K01_F01_E01"]) + ["k01_f00"]*len(config["K01_F00_E01"])):
            tree_name = f"{protein}_{ftype}_e01"
            if protein not in protein_list:
                protein_list.append(protein)
            filter_list.append(tree_name)
            protein_to_filter[protein] = tree_name
        
        # Add e00 trees
        for protein, ftype in zip(config["K01_F01_E00"] + config["K01_F00_E00"],
                                   ["k01_f01"]*len(config["K01_F01_E00"]) + ["k01_f00"]*len(config["K01_F00_E00"])):
            tree_name = f"{protein}_{ftype}_e00"
            if protein not in protein_list:
                protein_list.append(protein)
            filter_list.append(tree_name)
            protein_to_filter[protein] = tree_name
        
        # Sort lists
        protein_list.sort()
        filter_list.sort()
        
        # Write the config file
        with open(output[0], 'w') as f:
            f.write("# Configuration file for CD-HIT clustering (step 04)\n")
            f.write("# Auto-generated based on results from step 03 (e01/e00 filtering)\n\n")
            f.write("# List of all tree names to process\n")
            f.write("# Format: {PROTEIN}_{kegg_domain_filter}_{euk_filter}\n")
            f.write("# Examples: ADAT2_k01_f01_e01, CmoB_k01_f00_e01, ALKB1_k01_f01_e00, ArcS_k01_f00_e00\n")
            
            f.write("PROTEIN:\n")
            for p in protein_list:
                f.write(f"  - {p}\n")
            
            f.write("\nFILTER:\n")
            for tree in filter_list:
                f.write(f"  - {tree}\n")
            
            f.write("\n\nMAP: # Mapping from PROTEIN to TREE\n")
            for p in protein_list:
                if p in protein_to_filter:
                    f.write(f"  {p}: {protein_to_filter[p]}\n")
            
            f.write("\n\n# IQTree ultrafast bootstrap parameter\n")
            f.write("uf: 1000\n")

