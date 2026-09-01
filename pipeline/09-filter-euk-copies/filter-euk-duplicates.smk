# Copyright © 2024, Begüm Serra Büyüktarakçı
# All rights reserved

configfile: "filter-main.yaml"

results = "../../output/"
# Filter out TusC from the protein list
PROTEIN = [p for p in config["PROTEIN"] if not p.startswith("TusC_")]
FILTER_TO_PROTEIN = config["FILTER_TO_PROTEIN"]

# Create reverse mapping: full name -> base name
PROTEIN_TO_BASE = {v: k for k, v in FILTER_TO_PROTEIN.items()}

# Combined tree proteins (with different paths)
import glob
combined_files = glob.glob(f"{results}seqs/combined_tree/*_m01.fa")
COMBINED_PROTEINS = [f.split('/')[-1].replace('.fa', '') for f in combined_files]

wildcard_constraints:
    protein="|".join(PROTEIN),
    combined="|".join(COMBINED_PROTEINS)

rule all:
    input:
        expand(results + "seqs/homologs/{protein}_r00.fa", protein=PROTEIN),
        expand(results + "seqs/homologs/{protein}_r01.fa", protein=PROTEIN),
        expand(results + "seqs/combined_tree/{combined}_r00.fa", combined=COMBINED_PROTEINS),
        expand(results + "seqs/combined_tree/{combined}_r01.fa", combined=COMBINED_PROTEINS),
        "../11-tree-Cmodels/tree.yaml"
    run:
        import os
        # Cleanup: remove empty files after workflow completes
        for protein in PROTEIN:
            r00 = f"{results}seqs/homologs/{protein}_r00.fa"
            r01 = f"{results}seqs/homologs/{protein}_r01.fa"
            
            # Remove r00 if empty
            if os.path.exists(r00) and os.path.getsize(r00) == 0:
                os.remove(r00)
                print(f"Cleaned up empty {r00}")
            
            # Remove r01 if empty  
            if os.path.exists(r01) and os.path.getsize(r01) == 0:
                os.remove(r01)
                print(f"Cleaned up empty {r01}")
        
        # Cleanup combined tree files
        for protein in COMBINED_PROTEINS:
            r00 = f"{results}seqs/combined_tree/{protein}_r00.fa"
            r01 = f"{results}seqs/combined_tree/{protein}_r01.fa"
            
            if os.path.exists(r00) and os.path.getsize(r00) == 0:
                os.remove(r00)
                print(f"Cleaned up empty {r00}")
            
            if os.path.exists(r01) and os.path.getsize(r01) == 0:
                os.remove(r01)
                print(f"Cleaned up empty {r01}")

rule identify_euk_duplicates:
    """
    Identify eukaryotic duplicate sequences from the same species in the tree.
    Keep only the longest sequence and create a list of IDs to remove.
    """
    input:
        tree = results + "iqtree/{protein}.treefile",
        fasta = results + "seqs/homologs/{protein}.fa"
    output:
        remove_ids = results + "seqs/filtering/ids/r01/.done_{protein}.txt"
    params:
        base_protein = lambda wildcards: PROTEIN_TO_BASE[wildcards.protein],
        final_output = lambda wildcards: results + f"seqs/filtering/ids/r01/{PROTEIN_TO_BASE[wildcards.protein]}.txt"
    shell:
        """
        python final-clean-tree-euk.py {input.tree} {input.fasta} {params.final_output}
        touch {output.remove_ids}
        """

rule filter_fasta:
    """
    Remove duplicate eukaryotic sequences from fasta file.
    If IDs are removed, create _r01 suffix, otherwise _r00.
    Only the file with data will have content; the other will be removed.
    """
    input:
        fasta = results + "seqs/homologs/{protein}.fa",
        remove_ids = results + "seqs/filtering/ids/r01/.done_{protein}.txt"
    params:
        remove_file = lambda wildcards: results + f"seqs/filtering/ids/r01/{PROTEIN_TO_BASE[wildcards.protein]}.txt"
    output:
        fasta_r00 = results + "seqs/homologs/{protein}_r00.fa",
        fasta_r01 = results + "seqs/homologs/{protein}_r01.fa"
    run:
        from Bio import SeqIO
        import os
        
        # Read IDs to remove from the actual file (params.remove_file)
        remove_ids = set()
        if os.path.exists(params.remove_file):
            with open(params.remove_file) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        remove_ids.add(line)
        
        # Determine which file gets the data
        if len(remove_ids) > 0:
            # Sequences removed - r01 gets data, r00 is empty
            print(f"{wildcards.protein}: Removing {len(remove_ids)} sequences -> _r01 has data")
            data_file = str(output.fasta_r01)
            empty_file = str(output.fasta_r00)
        else:
            # No sequences removed - r00 gets data, r01 is empty
            print(f"{wildcards.protein}: No sequences to remove -> _r00 has data")
            data_file = str(output.fasta_r00)
            empty_file = str(output.fasta_r01)
        
        # Write filtered sequences to the data file
        kept_count = 0
        with open(data_file, 'w') as out_f:
            for record in SeqIO.parse(input.fasta, "fasta"):
                if record.id not in remove_ids:
                    SeqIO.write(record, out_f, "fasta")
                    kept_count += 1
        
        print(f"  Wrote {kept_count} sequences to {data_file}")
        
        # Create empty file for the other output (Snakemake requires it)
        open(empty_file, 'w').close()
        
        # Remove the temporary .done marker file
        done_file = str(input.remove_ids)
        if os.path.exists(done_file):
            os.remove(done_file)

rule identify_euk_duplicates_combined:
    """
    Identify eukaryotic duplicate sequences for combined tree proteins.
    """
    input:
        tree = results + "iqtree/combined_tree/{combined}.treefile",
        fasta = results + "seqs/combined_tree/{combined}.fa"
    output:
        remove_ids = results + "seqs/filtering/ids/r01/.done_combined_{combined}.txt"
    params:
        final_output = lambda wildcards: results + f"seqs/filtering/ids/r01/combined_{wildcards.combined}.txt"
    shell:
        """
        python final-clean-tree-euk.py {input.tree} {input.fasta} {params.final_output}
        touch {output.remove_ids}
        """

rule filter_fasta_combined:
    """
    Remove duplicate eukaryotic sequences from combined tree fasta files.
    """
    input:
        fasta = results + "seqs/combined_tree/{combined}.fa",
        remove_ids = results + "seqs/filtering/ids/r01/.done_combined_{combined}.txt"
    params:
        remove_file = lambda wildcards: results + f"seqs/filtering/ids/r01/combined_{wildcards.combined}.txt"
    output:
        fasta_r00 = results + "seqs/combined_tree/{combined}_r00.fa",
        fasta_r01 = results + "seqs/combined_tree/{combined}_r01.fa"
    run:
        from Bio import SeqIO
        import os
        
        # Read IDs to remove
        remove_ids = set()
        if os.path.exists(params.remove_file):
            with open(params.remove_file) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        remove_ids.add(line)
        
        # Determine which file gets the data
        if len(remove_ids) > 0:
            print(f"{wildcards.combined}: Removing {len(remove_ids)} sequences -> _r01 has data")
            data_file = str(output.fasta_r01)
            empty_file = str(output.fasta_r00)
        else:
            print(f"{wildcards.combined}: No sequences to remove -> _r00 has data")
            data_file = str(output.fasta_r00)
            empty_file = str(output.fasta_r01)
        
        # Write filtered sequences
        kept_count = 0
        with open(data_file, 'w') as out_f:
            for record in SeqIO.parse(input.fasta, "fasta"):
                if record.id not in remove_ids:
                    SeqIO.write(record, out_f, "fasta")
                    kept_count += 1
        
        print(f"  Wrote {kept_count} sequences to {data_file}")
        
        # Create empty file for the other output
        open(empty_file, 'w').close()
        
        # Remove the temporary .done marker file
        done_file = str(input.remove_ids)
        if os.path.exists(done_file):
            os.remove(done_file)

rule create_tree_config:
    """
    Create tree.yaml config file for the next step with r00/r01 suffixes.
    Includes both regular proteins and combined tree proteins.
    """
    input:
        expand(results + "seqs/homologs/{protein}_r00.fa", protein=PROTEIN),
        expand(results + "seqs/homologs/{protein}_r01.fa", protein=PROTEIN),
        expand(results + "seqs/combined_tree/{combined}_r00.fa", combined=COMBINED_PROTEINS),
        expand(results + "seqs/combined_tree/{combined}_r01.fa", combined=COMBINED_PROTEINS)
    output:
        "../11-tree-Cmodels/tree.yaml"
    run:
        import os
        
        # Determine which files actually exist (non-empty)
        protein_list = []
        filter_to_protein = {}
        
        # Process regular proteins
        for protein in PROTEIN:
            base_protein = PROTEIN_TO_BASE[protein]
            r00 = f"{results}seqs/homologs/{protein}_r00.fa"
            r01 = f"{results}seqs/homologs/{protein}_r01.fa"
            
            # Check which file has data
            if os.path.exists(r00) and os.path.getsize(r00) > 0:
                final_name = f"{protein}_r00"
                protein_list.append(final_name)
                filter_to_protein[base_protein] = final_name
            elif os.path.exists(r01) and os.path.getsize(r01) > 0:
                final_name = f"{protein}_r01"
                protein_list.append(final_name)
                filter_to_protein[base_protein] = final_name
        
        # Process combined proteins
        for protein in COMBINED_PROTEINS:
            # Use the combined protein name as the base (e.g., "MnmA-MTU1_m01")
            base_protein = protein
            r00 = f"{results}seqs/combined_tree/{protein}_r00.fa"
            r01 = f"{results}seqs/combined_tree/{protein}_r01.fa"
            
            # Check which file has data
            if os.path.exists(r00) and os.path.getsize(r00) > 0:
                final_name = f"{protein}_r00"
                protein_list.append(final_name)
                filter_to_protein[base_protein] = final_name
            elif os.path.exists(r01) and os.path.getsize(r01) > 0:
                final_name = f"{protein}_r01"
                protein_list.append(final_name)
                filter_to_protein[base_protein] = final_name
        
        # Write the config file
        with open(output[0], 'w') as f:
            f.write("# Configuration for phylogenetic tree construction\n")
            f.write("# Lists all proteins with eukaryotic duplicates filtered (r00 or r01)\n")
            f.write("# Includes both regular proteins and combined tree proteins\n")
            f.write("PROTEIN:\n")
            f.write("  [\"" + "\",\n  \"".join(protein_list) + "\"]\n\n")
            f.write("# Mapping from protein name to filter string (for file lookup)\n")
            f.write("FILTER_TO_PROTEIN:\n")
            for base, full in sorted(filter_to_protein.items()):
                f.write(f"  {base}: {full}\n")
        
        print(f"Created tree config with {len(protein_list)} proteins")
        print(f"  Regular proteins: {len(PROTEIN)}")
        print(f"  Combined proteins: {len(COMBINED_PROTEINS)}")
        print(f"Output: {output[0]}")


