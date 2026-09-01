import re
import sys

def extract_id_species(header):
    # Extract ID (first word after '>')
    id_ = header[1:].strip().split()[0]
    # If ID starts with 'EP', always extract species from ID, ignore brackets
    if id_.startswith('EP'):
        id_parts = id_.split('_')
        # Remove first part (EP...)
        id_parts = id_parts[1:]
        # Remove last part if it starts with 'P' and is alphanumeric
        if len(id_parts) > 1 and re.match(r'^P\w+$', id_parts[-1]):
            id_parts = id_parts[:-1]
        species = ' '.join(id_parts)
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

def fasta_to_csv(fasta_path, output_path):
    with open(fasta_path) as fasta, open(output_path, "w") as out_color:
        out_color.write("LABELS\n")
        out_color.write("SEPARATOR COMMA\n")
        out_color.write("\n")
        out_color.write("DATA\n")
        for line in fasta:
            if line.startswith(">"):
                id_, species = extract_id_species(line)
                # If id does not start with EP and species does not start with Ca, keep only first two words
                if not id_.startswith('EP') and not (species.startswith('Ca.') or species.startswith('Candidatus ')):
                    species = ' '.join(species.split()[:2])
                species = format_species_name(species)
                out_color.write(f"{id_},{species}\n")  # Always write ID, species may be blank

def main():
    if len(sys.argv) != 3:
        print("Usage: python fasta_to_csv.py input.fasta output.txt")
        sys.exit(1)
    input_fasta = sys.argv[1]
    output_txt = sys.argv[2]
    fasta_to_csv(input_fasta, output_txt)

if __name__ == "__main__":
    main()
