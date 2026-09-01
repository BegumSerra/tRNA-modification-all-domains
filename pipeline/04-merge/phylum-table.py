import pandas as pd
from Bio import SeqIO
import os
import sys
import re


def clean_id(record_id):
    if re.match(r'^EP\d+_', record_id) and re.search(r'_P\d+$', record_id):
        return True, re.sub(r'^EP\d+_|_P\d+$', '', record_id)
    else:
        return False, record_id


def extract_and_clean_ids(fasta_file):
    euk_species = []
    euk_ids = []
    ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        is_cleaned, cleaned_id = clean_id(record.id)
        if is_cleaned:
            euk_species.append(cleaned_id)
        else:
            ids.append(record.id)
    print("euk species", euk_species)
    print("ids", ids)

    for record in SeqIO.parse(fasta_file, "fasta"):
        if re.match(r'^EP\d+_', record.id) and re.search(r'_P\d+$', record.id):
            euk_ids.append(record.id)
    print("euk ids", euk_ids)

    return euk_species, euk_ids, ids


def get_protein_name_from_path(fasta_file):
    file_name = os.path.basename(fasta_file)
    protein_name = file_name.split('_')[0]
    return protein_name


def find_taxids(ids, id_to_taxid_file):
    id_to_taxid_df = pd.read_csv(id_to_taxid_file, sep='\t', header=None, names=['id', 'taxid'])
    id_to_taxid_dict = dict(zip(id_to_taxid_df['id'], id_to_taxid_df['taxid']))
    taxids = [id_to_taxid_dict.get(protein_id) for protein_id in ids if protein_id in id_to_taxid_dict]
    return taxids


def extract_taxonomy_from_string(taxonomy_str, rank_prefix):
    """Extract a specific rank from GTDB-formatted taxonomy string."""
    if pd.isna(taxonomy_str):
        return 'FALSE'
    parts = taxonomy_str.split(';')
    for part in parts:
        if part.startswith(rank_prefix):
            result = part.replace(rank_prefix, '').strip()
            return result if result else 'FALSE'
    return 'FALSE'


def find_taxonomy_info(taxids, taxid_to_species_file, euk_species):
    """Extract domain, phylum, order, family, and species from ncbi_taxonomy column."""
    taxid_to_species_df = pd.read_csv(taxid_to_species_file, sep='\t')
    
    # Extract taxonomic ranks from ncbi_taxonomy column
    taxid_to_species_df['domain'] = taxid_to_species_df['ncbi_taxonomy'].apply(
        lambda x: extract_taxonomy_from_string(x, 'd__'))
    taxid_to_species_df['phylum'] = taxid_to_species_df['ncbi_taxonomy'].apply(
        lambda x: extract_taxonomy_from_string(x, 'p__'))
    taxid_to_species_df['order'] = taxid_to_species_df['ncbi_taxonomy'].apply(
        lambda x: extract_taxonomy_from_string(x, 'o__'))
    taxid_to_species_df['family'] = taxid_to_species_df['ncbi_taxonomy'].apply(
        lambda x: extract_taxonomy_from_string(x, 'f__'))
    taxid_to_species_df['species'] = taxid_to_species_df['ncbi_taxonomy'].apply(
        lambda x: extract_taxonomy_from_string(x, 's__'))
    
    # If species is 'FALSE' or empty, use ncbi_organism_name as fallback
    def get_species_with_fallback(row):
        species = row['species']
        if species == 'FALSE' or pd.isna(species) or species == '':
            organism_name = row.get('ncbi_organism_name', 'FALSE')
            return organism_name if pd.notna(organism_name) and organism_name != '' else 'FALSE'
        return species
    
    taxid_to_species_df['species'] = taxid_to_species_df.apply(get_species_with_fallback, axis=1)
    
    # Create dictionaries for lookups
    species_taxid_dict = dict(zip(
        taxid_to_species_df['ncbi_species_taxid'], 
        zip(taxid_to_species_df['domain'], taxid_to_species_df['phylum'], 
            taxid_to_species_df['order'], taxid_to_species_df['family'], taxid_to_species_df['species'])
    ))
    taxid_dict = dict(zip(
        taxid_to_species_df['ncbi_taxid'], 
        zip(taxid_to_species_df['domain'], taxid_to_species_df['phylum'], 
            taxid_to_species_df['order'], taxid_to_species_df['family'], taxid_to_species_df['species'])
    ))

    taxonomy_info = []
    for taxid in taxids:
        if taxid in species_taxid_dict:
            taxonomy_info.append(species_taxid_dict[taxid])
        elif taxid in taxid_dict:
            taxonomy_info.append(taxid_dict[taxid])
        else:
            taxonomy_info.append(('FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE'))

    # Add eukaryote species (assuming they need domain and phylum too)
    for euk_sp in euk_species:
        taxonomy_info.append(('Eukaryota', 'FALSE', 'FALSE', 'FALSE', euk_sp if euk_sp else 'FALSE'))
    
    return taxonomy_info


def create_dataframe_from_taxonomy(taxonomy_info_list, ids, euk_ids):
    """Create a dataframe with domain, phylum, order, family, species, and protein IDs."""
    protein_name = get_protein_name_from_path(fasta_file)
    print(protein_name)

    all_ids = ids + euk_ids

    # Create dataframe from taxonomy info
    data = []
    for (domain, phylum, order, family, species), protein_id in zip(taxonomy_info_list, all_ids):
        data.append({
            'domain': domain,
            'phylum': phylum,
            'order': order,
            'family': family,
            'species': species,
            protein_name: protein_id
        })
    
    df = pd.DataFrame(data)
    return df

def transform_table(df):
    arcS_column = df.columns[-1]

    df_split = df.set_index(df.columns[:-1].tolist())[arcS_column].str.split(', ', expand=True).stack().reset_index(
        level=len(df.columns) - 1, drop=True).reset_index(name=arcS_column)

    new_columns = [arcS_column] + df.columns[:-1].tolist()
    df_split = df_split[new_columns]
    return df_split

def main(fasta_file, id_to_taxid_file, taxid_to_species_file, csv_file, output_csv_file):
    euk_species, euk_ids, ids = extract_and_clean_ids(fasta_file)
    taxids = find_taxids(ids, id_to_taxid_file)
    taxonomy_info = find_taxonomy_info(taxids, taxid_to_species_file, euk_species)
    df = create_dataframe_from_taxonomy(taxonomy_info, ids, euk_ids)

    transformed_df = transform_table(df)
    transformed_df.columns = ["acc", "domain", "phylum", "order", "family", "species"]
    transformed_df.to_csv(output_csv_file, sep='\t', index=False, header=True)
    return transformed_df

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    id_to_taxid_file = sys.argv[2]
    taxid_to_species_file = sys.argv[3]
    csv_file = sys.argv[4]
    output_csv_file = sys.argv[5]

    transformed_df = main(fasta_file, id_to_taxid_file, taxid_to_species_file, csv_file, output_csv_file)
    print(transformed_df)

