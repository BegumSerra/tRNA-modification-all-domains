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
        return '#FFA500'  # orange
    
    # Look up in the mapping file
    domain = domain_map.get(protein_id, '')
    if domain == 'Bacteria':
        return '#FFC0CB'  # pink
    elif domain == 'Archaea':
        return '#0000FF'  # blue
    elif domain == '':
        # If domain is missing, check genome_id
        genome_id = genome_map.get(protein_id, '')
        # GCA genomes from archaea dataset are archaea
        if genome_id.startswith('GCA'):
            return '#0000FF'  # blue (archaea)
        else:
            return '#808080'  # gray for unknown
    else:
        return '#808080'  # gray for unknown

def extract_id_species(header):
    # Extract ID (first word after '>')
    id_ = header[1:].strip().split()[0]
    # If ID starts with 'EP', always extract species from ID, ignore brackets
    if id_.startswith('EP'):
        id_parts = id_.split('_')
        # Remove first part (EP...)
        species_parts = id_parts[1:]
        # Remove last part if it starts with 'P' and is alphanumeric
        if len(species_parts) > 1 and re.match(r'^P\w+$', species_parts[-1]):
            species_parts = species_parts[:-1]
        # Join with underscores to keep original format
        species = '_'.join(species_parts)
    else:
        # Try to extract species name inside brackets
        match = re.search(r'\[([^\]]+)\]', header)
        if match:
            species = match.group(1)
        else:
            # Fallback: try to extract from ID as above
            id_parts = id_.split('_')
            if len(id_parts) >= 3:
                species = f"{id_parts[1]} {id_parts[2]}"
            elif len(id_parts) == 2:
                species = id_parts[1]
            else:
                species = ""
    return id_, species

def format_species_name(species):
    species = species.strip()
    if species.startswith("Candidatus "):
        # Replace 'Candidatus' with 'Ca.'
        rest = species[len("Candidatus "):].strip()
        parts = rest.split()
        if len(parts) > 1:
            # Ca. X. ...
            return f"Ca. {parts[0][0]}. {' '.join(parts[1:])}"
        elif len(parts) == 1:
            return f"Ca. {parts[0][0]}."
        else:
            return "Ca."
    elif species.startswith("Ca. "):
        # If already Ca., shorten the next word to initial and dot
        rest = species[len("Ca. "):].strip()
        parts = rest.split()
        if len(parts) > 1:
            return f"Ca. {parts[0][0]}. {' '.join(parts[1:])}"
        elif len(parts) == 1:
            return f"Ca. {parts[0][0]}."
        else:
            return "Ca."
    elif species:
        parts = species.split()
        if len(parts) > 1:
            # Shorten to first letter and dot, then the rest
            return f"{parts[0][0]}. {' '.join(parts[1:])}"
        else:
            return species
    else:
        return ""

def fasta_to_csv(fasta_path, output_path, mapping_file):
    # Load domain mapping
    domain_map, genome_map = load_domain_mapping(mapping_file)
    
    legend_block = """
LEGEND_TITLE Example
LEGEND_SCALE 1
LEGEND_POSITION_X 100
LEGEND_POSITION_Y 100
LEGEND_HORIZONTAL 1
LEGEND_SHAPES 1
LEGEND_COLORS rgba(255,0,0,0.5)
LEGEND_LABELS Example
LEGEND_SHAPE_SCALES 1
"""
    with open(fasta_path) as fasta, open(output_path, "w") as out:
        out.write("LABELS\n")
        out.write("SEPARATOR COMMA\n")
        out.write("DATA\n")
        for line in fasta:
            if line.startswith(">"):
                id_, species = extract_id_species(line)
                color = get_domain_color(id_, domain_map, genome_map)
                
                # For eukaryotes (EP prefix), species is in underscore format
                if id_.startswith('EP'):
                    # Extract the protein ID part (last segment starting with P)
                    id_parts = id_.split('_')
                    protein_id = id_parts[-1] if re.match(r'^P\w+$', id_parts[-1]) else ""
                    # Convert species from underscores to spaces
                    species_display = species.replace('_', ' ')
                    out.write(f"{id_},{color},{id_} {species_display}\n")
                else:
                    # For non-eukaryotes
                    # If species does not start with Ca, keep only first two words
                    if not (species.startswith('Ca.') or species.startswith('Candidatus ')):
                        species = ' '.join(species.split()[:2])
                    species_formatted = format_species_name(species)
                    # Write protein ID first, then long species name (unformatted)
                    out.write(f"{id_},{color},{id_} {species}\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python fasta_to_csv.py input.fasta mapping.tsv output.txt")
        sys.exit(1)
    input_fasta = sys.argv[1]
    mapping_file = sys.argv[2]
    output_txt = sys.argv[3]
    fasta_to_csv(input_fasta, output_txt, mapping_file)

if __name__ == "__main__":
    main()
