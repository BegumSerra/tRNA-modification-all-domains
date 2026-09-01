import sys
import re
import os

# Color palette for proteins - different colors for each protein
PROTEIN_COLORS = [
    "#E69F00",  # Orange
    "#56B4E9",  # Sky blue
    "#009E73",  # Bluish green
    "#F0E442",  # Yellow
    "#0072B2",  # Blue
    "#D55E00",  # Vermillion
    "#CC79A7",  # Reddish purple
    "#999999",  # Gray
]

def extract_protein_name(fasta_file):
    """Extract protein name from filename (e.g., MnmA_k01_f01... -> MnmA)"""
    filename = os.path.basename(fasta_file)
    # Get the part before the first underscore
    match = re.match(r'^([A-Za-z0-9]+)_', filename)
    if match:
        return match.group(1)
    return os.path.splitext(filename)[0]

def parse_fasta_ids(fasta_path):
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def main(fasta_files, output_file):
    # First pass: identify proteins and assign colors
    protein_color_map = {}
    protein_ids = {}  # Store all IDs for each protein
    color_idx = 0
    
    for fasta_file in fasta_files:
        protein_name = extract_protein_name(fasta_file)
        if protein_name not in protein_color_map:
            protein_color_map[protein_name] = PROTEIN_COLORS[color_idx % len(PROTEIN_COLORS)]
            color_idx += 1
        ids = parse_fasta_ids(fasta_file)
        protein_ids[protein_name] = ids
    
    # Create range entries - use first and last ID to define range for label placement
    # Format for range: ID1|ID2,label,position,color,style,size_factor,rotation
    all_entries = []
    for protein_name, ids in protein_ids.items():
        if len(ids) >= 2:
            # Use first and last ID to create a range - label will be in the middle
            first_id = ids[0]
            last_id = ids[-1]
            color = protein_color_map[protein_name]
            # position: -1 = external label
            all_entries.append(f"{first_id}|{last_id},{protein_name},-1,{color},bold,1.5,0\n")
        elif len(ids) == 1:
            # Single sequence - just label it directly
            color = protein_color_map[protein_name]
            all_entries.append(f"{ids[0]},{protein_name},-1,{color},bold,1.5,0\n")
    
    # Write output
    with open(output_file, "w") as out:
        out.write("DATASET_TEXT\n")
        out.write("SEPARATOR COMMA\n")
        out.write("DATASET_LABEL,Protein_Names\n")
        out.write("COLOR,#000000\n")
        out.write("DATA\n")
        for entry in all_entries:
            out.write(entry)
    
    print(f"Generated {output_file} with {len(all_entries)} range labels")
    print(f"Proteins: {', '.join(protein_color_map.keys())}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python dataset-text-protein-layer.py protein1.fasta protein2.fasta ... output_file")
        sys.exit(1)
    fasta_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    main(fasta_files, output_file)
