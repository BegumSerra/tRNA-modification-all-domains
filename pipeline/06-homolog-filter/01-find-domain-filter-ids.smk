# All rights reserved
configfile: "01-config.yaml"

# Paths
results = "../../output/"
data = "../../input/"
database = "../../database/"

# Load the protein list from the config file
PROTEIN = config["PROTEIN"]  # List of proteins

rule all:
    input:
        expand(results + "seqs/filtering/ids/f01/{protein}_wo_domain.txt", protein=PROTEIN),
        expand(results + "seqs/filtering/ids/f01/{protein}_w_extra_domain.txt", protein=PROTEIN),
        expand(results + "seqs/filtering/ids/f01/{protein}_merged.txt", protein=PROTEIN),
        "02-config.yaml"

rule filter_wo_domain:
    input:
        fasta=results + "seqs/filtering/fasta/{protein}_k01.fa",
        phylum=results + "table/phylum_table/{protein}_phylum.table.tsv",
        domain=results + "domain/blast/merged/{protein}.tsv",
        kegg_filter=results + "seqs/filtering/ids/k01/{protein}_non_kegg.txt"
    output:
        results + "seqs/filtering/ids/f01/{protein}_wo_domain.txt"
    shell:
        """
        python filter-copies-wo-domains.py {input.fasta} {input.phylum} {input.domain} {input.kegg_filter} {output}
        """

rule filter_w_extra_domain:
    input:
        fasta=results + "seqs/filtering/fasta/{protein}_k01.fa",
	    phylum=results + "table/phylum_table/{protein}_phylum.table.tsv",
        domain=results + "domain/blast/merged/{protein}.tsv",
        kegg_filter=results + "seqs/filtering/ids/k01/{protein}_non_kegg.txt"
    output:
        results + "seqs/filtering/ids/f01/{protein}_w_extra_domain.txt"
    shell:
        """
        python filter-copies-w-extra-domains.py {input.fasta} {input.phylum} {input.domain} {input.kegg_filter} {output}
        """

rule merge_ids:
    input:
        wo_domain = results + "seqs/filtering/ids/f01/{protein}_wo_domain.txt",
        w_extra_domain = results + "seqs/filtering/ids/f01/{protein}_w_extra_domain.txt",
        kegg_filter = results + "seqs/filtering/ids/k01/{protein}_non_kegg.txt"
    output:
        results + "seqs/filtering/ids/f01/{protein}_merged.txt"
    shell:
        """
        python merge-and-exclude-kegg.py {input.wo_domain} {input.w_extra_domain} {input.kegg_filter} {output}
        """

rule collect_skipped_proteins:
    input:
        merged_ids = expand(results + "seqs/filtering/ids/f01/{protein}_merged.txt", protein=PROTEIN)
    output:
        "02-config.yaml"
    run:
        import os
        import yaml

        # Proteins that have sequences to remove (non-empty merged files)
        proteins_with_f01 = []
        for merged_file in input.merged_ids:
            if os.path.getsize(merged_file) > 0:  # Check if file is not empty
                protein_name = os.path.basename(merged_file).replace('_merged.txt','')
                proteins_with_f01.append(protein_name)

        # Proteins that don't have sequences to remove (empty merged files or not in list)
        proteins_no_f01 = []
        for protein in PROTEIN:
            merged_file = results + f"seqs/filtering/ids/f01/{protein}_merged.txt"
            if os.path.exists(merged_file) and os.path.getsize(merged_file) == 0:
                proteins_no_f01.append(protein)

        # Save both lists in the config
        config_dict = {
            "F01": proteins_with_f01,
            "F00": proteins_no_f01
        }
        with open(output[0],"w") as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)