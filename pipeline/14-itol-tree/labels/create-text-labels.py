#!/usr/bin/env python3
"""
Create iTOL text dataset with domain-based colors.
Format: ID,label,position,color,style,size_factor,rotation
"""

import re
import sys
import pandas as pd

def load_domain_mapping(mapping_file):
    """Load protein_id to domain and genome_id mapping from TSV file."""
    df = pd.read_csv(mapping_file, sep='\t')
    # Create dictionaries: protein_id -> domain and protein_id -> genome_id
    domain_map = dict(zip(df['protein_id'], df['domain']))
    genome_map = dict(zip(df['protein_id'], df['genome_id']))
    return domain_map, genome_map
  
def get_domain_color(protein_id, domain_map, genome_map):
    """Return color based on domain: Eukaryota=orange, Bacteria=pink, Archaea=blue."""
    # Check if it's a eukaryote based on EP prefix
    if protein_id.startswith('EP'):
        return '#FE6100'  # orange
    
    # Look up in the mapping file
    domain = domain_map.get(protein_id, '')
    if domain == 'Bacteria':
        return '#DC267F' # pink
    elif domain == 'Archaea':
        return '#648FFF'  # blue
    elif domain == '':
        # If domain is missing, check genome_id
        genome_id = genome_map.get(protein_id, '')
        # GCA genomes from archaea dataset are archaea
        if genome_id.startswith('GCA'):
            return '#648FFF'  # blue (archaea)
        else:
            return '#808080'  # gray for unknown
    else:
        return '#808080'  # gray for unknown

def extract_species_label(header, protein_id):
    """Extract species name for the label with protein ID."""
    # For eukaryotes (EP prefix)
    if protein_id.startswith('EP'):
        id_parts = protein_id.split('_')
        # Get the protein ID part (last segment starting with P)
        prot_id = id_parts[-1] if re.match(r'^P\w+$', id_parts[-1]) else ""
        # Get species parts (between EP... and P...)
        species_parts = id_parts[1:]
        # Remove last part if it starts with 'P' and is alphanumeric
        if len(species_parts) > 1 and re.match(r'^P\w+$', species_parts[-1]):
            species_parts = species_parts[:-1]
        # Convert underscores to spaces
        species = ' '.join(species_parts)
        return f"{prot_id} {species}"
    else:
        # For bacteria/archaea, extract from brackets
        match = re.search(r'\[([^\]]+)\]', header)
        if match:
            species = match.group(1)
            # Keep only first two words if not Candidatus
            if not species.startswith('Candidatus ') and not species.startswith('Ca.'):
                species = ' '.join(species.split()[:2])
            return f"{protein_id} {species}"
        else:
            # Fallback: try to extract from ID
            id_parts = protein_id.split('_')
            if len(id_parts) >= 3:
                return f"{protein_id} {id_parts[1]} {id_parts[2]}"
            elif len(id_parts) == 2:
                return f"{protein_id} {id_parts[1]}"
            else:
                return protein_id

def fasta_to_itol_text(fasta_path, mapping_file, output_path):
    """Convert FASTA to iTOL text dataset format."""
    # Load domain mapping
    domain_map, genome_map = load_domain_mapping(mapping_file)
    
    with open(fasta_path) as fasta, open(output_path, "w") as out:
        # Write header
        out.write("DATASET_TEXT\n")
        out.write("SEPARATOR COMMA\n")
        out.write("DATASET_LABEL,Species labels\n")
        out.write("COLOR,#000000\n\n")
        
        # Optional settings
        out.write("MARGIN,0\n")
        out.write("SHOW_INTERNAL,0\n")
        out.write("LABELS_BELOW,1\n")
        out.write("SIZE_FACTOR,1\n\n")
        
        # Data section
        out.write("DATA\n")
        out.write("#ID,label,position,color,style,size_factor,rotation\n")
        
        for line in fasta:
            if line.startswith(">"):
                # Extract protein ID (first word after >)
                protein_id = line[1:].strip().split()[0]
                
                # Get domain color
                color = get_domain_color(protein_id, domain_map, genome_map)
                
                # Extract species label
                species_label = extract_species_label(line, protein_id)
                
                # Format: ID,label,position,color,style,size_factor,rotation
                # position=-1 for external label, style=normal, size=1, rotation=0
                out.write(f"{protein_id},{species_label},-1,{color},normal,1,0\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python create-text-labels.py input.fasta mapping.tsv output.txt")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    mapping_file = sys.argv[2]
    output_txt = sys.argv[3]
    
    fasta_to_itol_text(input_fasta, mapping_file, output_txt)
    print(f"Created iTOL text dataset: {output_txt}")

if __name__ == "__main__":
    main()
