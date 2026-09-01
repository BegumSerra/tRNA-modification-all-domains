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
        # Read the domain file
        domain_data = pd.read_csv(domain_file, sep="\t", header=None)
        domain_data = domain_data[[0, 4]]  # Keep only the 1st and 5th columns
        domain_data.columns = ['Protein_ID', 'Domain']
        
        # Filter to only PFAM domains (starting with 'PF')
        domain_data = domain_data[domain_data['Domain'].str.startswith('PF', na=False)]

        # Group by Protein_ID and combine domains into a single set per protein
        grouped_domains = domain_data.groupby('Protein_ID')['Domain'].apply(
            lambda x: ','.join(x.dropna().unique())
        ).reset_index()
        grouped_domains.columns = ['Protein_ID', 'Domains']
        return grouped_domains
    except Exception as e:
        print("Error loading domain file: {}".format(e))
        sys.exit(1)


def parse_domains(domains_str):
    # Split domains and standardize domain IDs
    domains = [d.strip().upper() for d in domains_str.split(',')]
    print(domains)
    return domains

def get_domain_counts(domains_list):
    # Count the occurrences of each domain
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
    domain_data = load_domain_info(domain_file)
    kegg_filtered_ids = load_kegg_filter_ids(kegg_filter_file)
    
    # Exclude proteins already filtered by KEGG
    protein_ids = [pid for pid in protein_ids if pid not in kegg_filtered_ids]
    print("Proteins after KEGG filtering: {}".format(len(protein_ids)))

    # Filter domain data and phylum table to only proteins in the tree
    domain_data = domain_data[domain_data['Protein_ID'].isin(protein_ids)]
    phylum_table = phylum_table[phylum_table['acc'].isin(protein_ids)]

    # Remove proteins with missing or empty domain information
    domain_data = domain_data.dropna(subset=['Domains'])
    domain_data = domain_data[domain_data['Domains'].str.strip().astype(bool)]  # Remove empty strings

    # Parse domains and get domain counts for each protein
    domain_data['Domains_List'] = domain_data['Domains'].apply(parse_domains)
    domain_data['Domain_Counts'] = domain_data['Domains_List'].apply(get_domain_counts)

    # Determine the regular domain profile (with counts)
    # Count the frequency of each domain profile
    domain_profile_counts = domain_data['Domain_Counts'].value_counts()
    if domain_profile_counts.empty:
        print("No domain profiles found. Exiting.")
        sys.exit(1)

    regular_domain_profile = domain_profile_counts.index[0]
    print("\nRegular Domain Profile (Most Common Domain Counts):", regular_domain_profile)

    # Group proteins by species
    species_dict = defaultdict(list)
    for _, row in phylum_table.iterrows():
        species_dict[row['species']].append(row['acc'])

    # For species with multiple copies
    to_remove = []
    for species, proteins in species_dict.items():
        if len(proteins) > 1:
            # Get domain info for these proteins
            species_domain_subset = domain_data[domain_data['Protein_ID'].isin(proteins)]
            if species_domain_subset.empty:
                continue

            # Check if any protein in this species has the regular domain profile
            has_regular_domain = species_domain_subset['Domain_Counts'].isin([regular_domain_profile]).any()
            if not has_regular_domain:
                # No copy with the regular domain profile, do not remove any proteins
                #print("\nSpecies '{{species}}' does not have any protein with the regular domain profile. Skipping.")
                continue
            else:
                print("\nSpecies '{{species}}' has at least one protein with the regular domain profile.")

            # Identify proteins to remove based on the specified logic
            for _, row in species_domain_subset.iterrows():
                protein_id = row['Protein_ID']
                protein_domain_counts = row['Domain_Counts']

                # Compare domain counts
                if protein_domain_counts == regular_domain_profile:
                    # Domains match the regular domain profile
                    continue  # Keep the protein

                # Condition 1: Has completely different domains (no overlap with regular profile)
                # Remove if protein has domains but NONE of them are in the regular domain profile
                if protein_domain_counts and not any(domain in regular_domain_profile for domain in protein_domain_counts):
                    # Only remove if this species has another copy with regular domains
                    if has_regular_domain:
                        to_remove.append(protein_id)
                        print(f"Protein '{protein_id}' has completely different domains (no overlap) and species has regular copy. Marked for removal.")
                    continue

                # Condition 2: Has extra RARE domains (domains not in the regular profile)
                # Only flag if there are extra domains beyond the regular ones
                extra_domains = protein_domain_counts - regular_domain_profile
                if len(extra_domains) > 0:
                    # This protein has all regular domains plus additional rare domains
                    if has_regular_domain:
                        to_remove.append(protein_id)
                        print(f"Protein '{protein_id}' has extra rare domains: {list(extra_domains.keys())}. Marked for removal.")
                    continue

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

