# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

import os

configfile: "04-config.yaml"

results = "../../output/"
data = "../../input/"
database = "../../database/"

# Load configuration
PROTEIN = config["PROTEIN"]
FILTER = config["FILTER"]
protein_to_filter = config["MAP"]
filter_to_protein = {filter_val: protein_key for protein_key, filter_val in protein_to_filter.items()}

wildcard_constraints:
    filter = "|".join(FILTER)

def get_protein_id_from_filter(filter):
    if filter not in filter_to_protein:
        raise ValueError(f"Filter '{filter}' not found in mapping.")
    return filter_to_protein[filter]

def get_ref_seq_path(protein, domain):
    path = f"{data}ref_seqs/{domain}/{protein}.fa"
    if os.path.exists(path):
        return path
    return f"{data}ref_seqs/{domain}/.empty_{protein}.fa"

rule all:
    input:
        expand(results + "seqs/filtering/cdhit/{filter}.clstr", filter=FILTER),
        expand(results + "seqs/filtering/fasta/{filter}_c01.fa", filter=FILTER),
        expand(results + "seqs/filtering/fasta/{filter}_c00.fa", filter=FILTER),
        "05-summary.yaml",
        results + "table/filtering_summary.tsv"

rule cdhit:
    input:
        results + "seqs/filtering/fasta/{filter}.fa"
    output:
        results + "seqs/filtering/cdhit/{filter}.cdhit.fa",
        results + "seqs/filtering/cdhit/{filter}.clstr"
    params:
        prefix=results + "seqs/filtering/cdhit/{filter}"
    shell:
        """
        cd-hit -d 0 -i {input} -o {params.prefix} -c 0.95
        mv {params.prefix} {output[0]}
        """

rule add_back_seqs:
    input:
        original=results + "seqs/filtering/fasta/{filter}.fa",
        cdhit=results + "seqs/filtering/cdhit/{filter}.cdhit.fa",
        clstr=results + "seqs/filtering/cdhit/{filter}.clstr"
    output:
        c01=results + "seqs/filtering/fasta/{filter}_c01.fa",
        c00=results + "seqs/filtering/fasta/{filter}_c00.fa"
    run:
        from Bio import SeqIO
        import os
        
        # Count sequences before clustering
        original_count = sum(1 for _ in SeqIO.parse(input.original, "fasta"))
        
        # Run add-back script
        shell("python add-back-cdhit.py {input.original} {input.cdhit} {input.clstr} {output.c01}")
        
        # Count after add-back
        clustered_count = sum(1 for _ in SeqIO.parse(output.c01, "fasta"))
        
        # Determine if clustering happened
        if clustered_count < original_count:
            # Clustering reduced sequences - keep c01, make c00 empty marker
            shell("touch {output.c00}")
        else:
            # No clustering - copy to c00, make c01 empty marker
            shell("cp {output.c01} {output.c00}")
            shell("echo '# No clustering occurred' > {output.c01}")

rule collect_clustering_results:
    input:
        c01_files = expand(results + "seqs/filtering/fasta/{filter}_c01.fa", filter=FILTER),
        c00_files = expand(results + "seqs/filtering/fasta/{filter}_c00.fa", filter=FILTER)
    output:
        "05-summary.yaml"
    run:
        import yaml
        import os
        
        # Identify which filters have c01 (clustering) vs c00 (no clustering)
        c01_filters = []
        c00_filters = []
        
        for filter_name in FILTER:
            c01_file = results + f"seqs/filtering/fasta/{filter_name}_c01.fa"
            c00_file = results + f"seqs/filtering/fasta/{filter_name}_c00.fa"
            
            # Check file sizes - non-empty file is the real one
            c01_size = os.path.getsize(c01_file) if os.path.exists(c01_file) else 0
            c00_size = os.path.getsize(c00_file) if os.path.exists(c00_file) else 0
            
            if c01_size > c00_size:
                c01_filters.append(filter_name)
                # Delete the empty c00 file
                if os.path.exists(c00_file) and c00_size < 100:
                    os.remove(c00_file)
            else:
                c00_filters.append(filter_name)
                # Delete the empty c01 file
                if os.path.exists(c01_file) and c01_size < 100:
                    os.remove(c01_file)
        
        # Create config
        config_dict = {
            "C01": c01_filters,
            "C00": c00_filters,
        }
        
        with open(output[0], "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

rule filtering_summary:
    input:
        "05-summary.yaml"
    output:
        results + "table/filtering_summary.tsv"
    params:
        results_fasta=results + "seqs/filtering/fasta",
        blast_merged=results + "seqs/blast/merged"
    shell:
        "python filtering_summary.py {params.results_fasta} {params.blast_merged} {output}"