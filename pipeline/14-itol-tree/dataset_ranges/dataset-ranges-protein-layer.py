import sys
import re
import os

# Color palette for proteins - will cycle through if more proteins
PROTEIN_COLORS = ["rgba(128,128,128,0.5)","rgba(128,128,128,1)"]

def extract_protein_name(fasta_file):
    """Extract protein name from filename (e.g., MnmA_k01_f01... -> MnmA)"""
    filename = os.path.basename(fasta_file)
    # Get the part before the first underscore
    match = re.match(r'^([A-Za-z0-9]+)_', filename)
    if match:
        return match.group(1)
    return os.path.splitext(filename)[0]

def generate_legend_block(protein_color_map):
    """Generate legend block based on proteins used"""
    proteins = list(protein_color_map.keys())
    colors = [protein_color_map[p] for p in proteins]
    n = len(proteins)
    
    shapes = "\t".join(["1"] * n)
    colors_str = "\t".join(colors)
    labels_str = "\t".join(proteins)
    scales = "\t".join(["1"] * n)
    
    return f"""LEGEND_TITLE\tProtein_Name
LEGEND_SHAPES\t{shapes}
LEGEND_COLORS\t{colors_str}
LEGEND_LABELS\t{labels_str}
"""

def parse_fasta_ids(fasta_path):
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def main(fasta_files, color_file):
    # First pass: identify proteins and assign colors
    protein_color_map = {}
    color_idx = 0
    
    for fasta_file in fasta_files:
        protein_name = extract_protein_name(fasta_file)
        if protein_name not in protein_color_map:
            protein_color_map[protein_name] = PROTEIN_COLORS[color_idx % len(PROTEIN_COLORS)]
            color_idx += 1
    
    # Generate legend
    legend_block = generate_legend_block(protein_color_map)
    
    # Second pass: collect entries
    all_entries = []
    for fasta_file in fasta_files:
        protein_name = extract_protein_name(fasta_file)
        color = protein_color_map[protein_name]
        ids = parse_fasta_ids(fasta_file)
        for id_ in ids:
            # Format: node_id\tstart_node\tcolor (DATASET_RANGE doesn't support labels)
            all_entries.append(f"{id_}\t{id_}\t{color}\n")
    
    # Write output
    with open(color_file, "w") as out_color:
        out_color.write("DATASET_RANGE\n")
        out_color.write("SEPARATOR TAB\n")
        out_color.write("DATASET_LABEL\tProtein_Layer\n")
        out_color.write("COLOR\t#ff0000\n")
        out_color.write(legend_block)
        out_color.write("DATA\n")
        for entry in all_entries:
            out_color.write(entry)
    
    print(f"Generated {color_file} with {len(all_entries)} entries")
    print(f"Proteins: {', '.join(protein_color_map.keys())}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python dataset-ranges-protein-layer.py protein1.fasta protein2.fasta ... output_file")
        sys.exit(1)
    fasta_files = sys.argv[1:-1]
    color_file = sys.argv[-1]
    main(fasta_files, color_file)
