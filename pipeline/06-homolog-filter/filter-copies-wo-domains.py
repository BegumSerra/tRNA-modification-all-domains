from Bio import SeqIO
import pandas as pd
from collections import defaultdict, Counter
import sys

def load_fasta(fasta_file):
    """Load FASTA file and extract protein IDs"""
    try:
        protein_ids = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_ids.append(record.id)
        return protein_ids
    except Exception as e:
        print("Error loading FASTA file: {{{0}}}".format(e))
        sys.exit(1)

def load_phylum_table(phylum_file):
    try:
        phylum_data = pd.read_csv(phylum_file, sep="\t")
        phylum_data = phylum_data[['acc', 'species']]
        return phylum_data
    except Exception as e:
        print("Error loading phylum table: {}".format(e))
        sys.exit(1)

def load_domain_info(domain_file):
    try:
        domain_info = pd.read_csv(domain_file, sep="\t", header=None)
        domain_info = domain_info[[0, 4]]  # Adjust columns if needed
        domain_info.columns = ['Protein_ID', 'Domains']
        
        # Filter to only PFAM domains (starting with 'PF')
        domain_info = domain_info[domain_info['Domains'].str.startswith('PF', na=False)]
        
        # Group by Protein_ID and combine PFAM domains
        domain_info = domain_info.groupby('Protein_ID')['Domains'].apply(
            lambda x: ','.join(x.dropna().unique())
        ).reset_index()
        
        return domain_info
    except Exception as e:
        print("Error loading domain file: {}".format(e))
        sys.exit(1)

def parse_domains(domains_str):
    if not domains_str.strip():
        return []
    domains = [d.strip().upper() for d in domains_str.split(',')]
    return domains

def get_domain_counts(domains_list):
    return Counter(domains_list)

def load_kegg_filter_ids(kegg_filter_file):
    """Load protein IDs that were filtered out by KEGG filtering"""
    try:
        with open(kegg_filter_file, 'r') as f:
            kegg_filtered_ids = set(line.strip() for line in f if line.strip())
        return kegg_filtered_ids
    except Exception as e:
        print("Error loading KEGG filter file: {}".format(e))
        sys.exit(1)

def main(fasta_file, phylum_file, domain_file, kegg_filter_file, output_removed_ids_file):
    # Load data
    protein_ids = load_fasta(fasta_file)
    print("Total proteins in FASTA: {}".format(len(protein_ids)))
    phylum_table = load_phylum_table(phylum_file)
    domain_info = load_domain_info(domain_file)
    kegg_filtered_ids = load_kegg_filter_ids(kegg_filter_file)
    
    # Exclude proteins already filtered by KEGG
    protein_ids = [pid for pid in protein_ids if pid not in kegg_filtered_ids]
    print("Proteins after KEGG filtering: {}".format(len(protein_ids)))

    # Create domain_data DataFrame with all proteins in the tree
    domain_data = pd.DataFrame({'Protein_ID': protein_ids})
    # Merge with domain_info to get Domains
    domain_data = domain_data.merge(domain_info, on='Protein_ID', how='left')
    # Replace NaN in 'Domains' with empty string
    domain_data['Domains'] = domain_data['Domains'].fillna('')

    # Parse domains and get domain counts for each protein
    domain_data['Domains_List'] = domain_data['Domains'].apply(parse_domains)
    domain_data['Domain_Counts'] = domain_data['Domains_List'].apply(get_domain_counts)

    # Determine the regular domain profile (with counts)
    # Exclude proteins with empty domain lists
    domain_data_non_empty = domain_data[domain_data['Domains_List'].map(bool)]
    domain_profile_counts = domain_data_non_empty['Domain_Counts'].value_counts()
    if domain_profile_counts.empty:
        print("No domain profiles found. Exiting.")
        sys.exit(1)

    regular_domain_profile = domain_profile_counts.index[0]
    print("\nRegular Domain Profile (Most Common Domain Counts):", regular_domain_profile)

    # Filter phylum table to include only proteins in the tree
    phylum_table = phylum_table[phylum_table['acc'].isin(protein_ids)]

    # Map proteins to species
    species_dict = defaultdict(list)
    for _, row in phylum_table.iterrows():
        species_dict[row['species']].append(row['acc'])

    # For species with multiple copies
    to_remove = []
    for species, proteins in species_dict.items():
        if len(proteins) > 1:
            # Get domain info for these proteins
            species_domain_subset = domain_data[domain_data['Protein_ID'].isin(proteins)]

            # Check if any protein in this species has the regular domain profile
            has_regular_domain = species_domain_subset['Domain_Counts'].isin([regular_domain_profile]).any()
            if not has_regular_domain:
                # No copy with the regular domain profile, do not remove any proteins
                print("\nSpecies '{{species}}' does not have any protein with the regular domain profile. Skipping.")
                continue
            else:
                print("\nSpecies '{{species}}' has at least one protein with the regular domain profile.")

            # Identify proteins to remove (those without domain information)
            for _, row in species_domain_subset.iterrows():
                protein_id = row['Protein_ID']
                domains_list = row['Domains_List']
                if not domains_list:
                    # Protein has no domain information
                    to_remove.append(protein_id)
                    print("Protein '{{protein_id}}' has no domain information. Marked for removal.")

    # Remove duplicates and sort
    unique_to_remove = sorted(set(to_remove))

    # Save to file
    with open(output_removed_ids_file, "w") as f:
        for protein_id in unique_to_remove:
            f.write(protein_id + "\n")
    print("\nUnique IDs of redundant sequences saved to {{output_removed_ids_file}}")

# Run the script
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script_name.py fasta_file phylum_file domain_file kegg_filter_file output_removed_ids_file")
        sys.exit(1)

    fasta_file = sys.argv[1]
    phylum_file = sys.argv[2]
    domain_file = sys.argv[3]
    kegg_filter_file = sys.argv[4]
    output_removed_ids_file = sys.argv[5]
    main(fasta_file, phylum_file, domain_file, kegg_filter_file, output_removed_ids_file)
