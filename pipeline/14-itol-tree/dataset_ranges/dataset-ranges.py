import random
import re
import sys
import csv

DOMAIN_COLORS = {
    'Bacteria': 'rgba(245,140,195,0.3)',
    'Eukarya': 'rgba(255,184,113,0.3)',
    'Archaea': 'rgba(100,149,237,0.3)',
}

PHYLA_COLORS = {
    'Cyano': 'rgba(76,142,115,0.6)',
    'Alpha': 'rgba(255,76,0,0.6)',
    'Chlamydia': 'rgba(247,255,71,0.6)',
}

CHLAMYDIA_GROUP = {
    'Chlamydia trachomatis',
    'Candidatus Protochlamydia amoebophila',
    'Waddlia chondrophila',
    'Simkania negevensis',
    'Chlamydiifrater phoenicopteri'
}

CYANO_GROUP = {
    'Nostoc punctiforme',
    'Planktothrix agardhii',
    'Prochlorococcus marinus',
    'Thermostichus lividus',
    'Gloeobacter violaceus'
}

ALPHA_GROUP = {
    'Neorickettsia sennetsu',
    'Gluconobacter oxydans',
    'Parvularcula bermudensis',
    'Phenylobacterium zucineum',
    'Zymomonas mobilis'
}

FIXED_COLORS = {
    "Thermoproteota": "rgba(89,75,124,0.3)", #594b7c
    "Euryarchaeota": "rgba(171,141,228,0.3)", #ab8de4
    "Candidatus Thermoplasmatota": "rgba(173,208,231,0.3)", #add0e7
    "Nitrososphaerota": "rgba(103,139,153,0.3)", #678b99
    "Asgard": "rgba(85,125,198,0.3)" #557dc6
}

ASGARD_GROUP = {
    "Candidatus Odinarchaeota",
    "Candidatus Lokiarchaeota",
    "Candidatus Heimdallarchaeota",
    "Candidatus Thorarchaeota",
}

# Map colors to legend labels (color -> label)
COLOR_TO_LABEL = {
    'rgba(245,140,195,0.3)': ('Bacteria', 'rgba(245,140,195,0.7)'),
    'rgba(255,184,113,0.3)': ('Eukarya', 'rgba(255,184,113,0.7)'),
    'rgba(100,149,237,0.3)': ('Archaea', 'rgba(100,149,237,0.7)'),
    'rgba(89,75,124,0.3)': ('Thermoproteota', 'rgba(89,75,124,0.7)'),
    'rgba(171,141,228,0.3)': ('Euryarchaeota', 'rgba(171,141,228,0.7)'),
    'rgba(173,208,231,0.3)': ('Thermoplasmatota', 'rgba(173,208,231,0.7)'),
    'rgba(103,139,153,0.3)': ('Nitrososphaerota', 'rgba(103,139,153,0.7)'),
    'rgba(85,125,198,0.3)': ('Asgard', 'rgba(85,125,198,0.7)'),
    'rgba(76,142,115,0.6)': ('Cyano', 'rgba(76,142,115,0.7)'),
    'rgba(255,76,0,0.6)': ('Alpha', 'rgba(255,76,0,0.6)'),
    'rgba(247,255,71,0.3)': ('Chlamydia', 'rgba(247,255,71,0.7)'),
}

def generate_legend_block(used_colors):
    """Generate legend block based on colors actually used in the data"""
    labels = []
    legend_colors = []
    
    for color in used_colors:
        if color in COLOR_TO_LABEL:
            label, legend_color = COLOR_TO_LABEL[color]
            if label not in labels:  # Avoid duplicates
                labels.append(label)
                legend_colors.append(legend_color)
    
    if not labels:
        return ""
    
    n = len(labels)
    shapes = " ".join(["1"] * n)
    colors_str = " ".join(legend_colors)
    labels_str = " ".join(labels)
    scales = " ".join(["1"] * n)
    
    return f"""
LEGEND_TITLE Taxonomy
LEGEND_SCALE 1
LEGEND_POSITION_X 100
LEGEND_POSITION_Y 100
LEGEND_HORIZONTAL 1
LEGEND_SHAPES {shapes}
LEGEND_COLORS {colors_str}
LEGEND_LABELS {labels_str}
LEGEND_SHAPE_SCALES {scales}
"""

def normalize_name(name: str) -> str:
    name = name.strip().lower()
    name = re.sub(r'\s+', ' ', name)
    return name

def sanitize_name(name: str) -> str:
    return re.sub(r'[^a-zA-Z0-9 _\-\.]', '', name)

def parse_fasta_ids(fasta_path):
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def parse_phylum_table(table_path):
    phylum_info = {}
    with open(table_path) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                seq_id = parts[0]
                domain = parts[1]
                phylum = parts[2]
                species = parts[5]
                phylum_info[seq_id] = (domain, phylum, species)
    return phylum_info

def assign_color(domain, phylum, species, protein=None):
    s_norm = normalize_name(species)
    phylum_norm = normalize_name(phylum)
    # For MnmG and MnmE, color Chlamydia as Bacteria
    if protein in {"MnmG", "MnmE"} and any(sub in s_norm for sub in [
        "chlamydia trachomatis", "protochlamydia amoebophila",
        "waddlia chondrophila", "simkania negevensis",
        "chlamydiifrater phoenicopteri"]):
        return DOMAIN_COLORS['Bacteria']
    # Color by species name for Cyano, Alpha, Chlamydia
    if s_norm in {normalize_name(s) for s in CYANO_GROUP}:
        return PHYLA_COLORS['Cyano']
    if s_norm in {normalize_name(s) for s in ALPHA_GROUP}:
        return PHYLA_COLORS['Alpha']
    if s_norm in {normalize_name(s) for s in CHLAMYDIA_GROUP}:
        return PHYLA_COLORS['Chlamydia']
    # For Asgard group, use Asgard color if phylum matches (normalized)
    if domain == 'Archaea':
        asgard_norm = {normalize_name(s) for s in ASGARD_GROUP}
        if phylum_norm in asgard_norm:
            return FIXED_COLORS["Asgard"]
        if phylum in FIXED_COLORS:
            return FIXED_COLORS[phylum]
    # If only archaea species are present, use FIXED_COLORS for phylum
    if domain == 'Archaea' and phylum in FIXED_COLORS:
        return FIXED_COLORS[phylum]
    # If not matched above, fallback to domain color
    if domain in DOMAIN_COLORS:
        return DOMAIN_COLORS[domain]
    # If not matched above, fallback to phylum color
    if phylum in PHYLA_COLORS:
        return PHYLA_COLORS[phylum]
    if phylum in FIXED_COLORS:
        return FIXED_COLORS[phylum]
    random.seed(phylum)
    # Fallback: use a random rgba color
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)
    return f"#000000 #000000 rgba({r},{g},{b},0.3)"

def main(fasta_file, phylum_table, color_file, protein=None, debug=False, mrca_ids_file=None):
    ids = parse_fasta_ids(fasta_file)
    phylum_info = parse_phylum_table(phylum_table)
    # Check if all domains are Archaea
    all_domains = set()
    for id_ in ids:
        if id_ in phylum_info:
            domain, phylum, species = phylum_info[id_]
            all_domains.add(domain)
    only_archaea = (len(all_domains) == 1 and 'Archaea' in all_domains)
    mrca_ids = set()
    if mrca_ids_file and protein in {"MnmG", "MnmE", "TilS"}:
        try:
            with open(mrca_ids_file) as f:
                for line in f:
                    mrca_ids.add(line.strip())
        except Exception:
            pass
    
    # First pass: collect all colors used
    id_colors = []
    used_colors = set()
    for id_ in ids:
        if id_ in phylum_info:
            domain, phylum, species = phylum_info[id_]
            species_clean = sanitize_name(species)
            s_norm = normalize_name(species)
            # Color by species name for Cyano, Alpha, Chlamydia
            if s_norm in {normalize_name(s) for s in CYANO_GROUP}:
                color = PHYLA_COLORS['Cyano']
            elif s_norm in {normalize_name(s) for s in ALPHA_GROUP}:
                color = PHYLA_COLORS['Alpha']
            elif s_norm in {normalize_name(s) for s in CHLAMYDIA_GROUP}:
                color = PHYLA_COLORS['Chlamydia']
            elif protein in {"MnmG", "MnmE", "TilS"} and id_ in mrca_ids:
                color = "rgba(76,142,115,0.3)"  # darker green, alpha 0.6
            else:
                if only_archaea:
                    phylum_norm = normalize_name(phylum)
                    asgard_norm = {normalize_name(s) for s in ASGARD_GROUP}
                    if phylum_norm in asgard_norm:
                        color = FIXED_COLORS["Asgard"]
                    elif phylum in FIXED_COLORS:
                        color = FIXED_COLORS[phylum]
                    else:
                        color = DOMAIN_COLORS['Archaea']
                else:
                    # Use PHYLA_COLORS if phylum matches, else DOMAIN_COLORS
                    if phylum in PHYLA_COLORS:
                        color = PHYLA_COLORS[phylum]
                    else:
                        # Normalize domain name for robust matching
                        domain_clean = domain.strip().lower()
                        if domain_clean in ['eukarya', 'eukaryota', 'eukaryotes']:
                            color = DOMAIN_COLORS['Eukarya']
                        elif domain_clean in ['bacteria', 'bacterial']:
                            color = DOMAIN_COLORS['Bacteria']
                        elif domain_clean in ['archaea', 'archaeal']:
                            color = DOMAIN_COLORS['Archaea']
                        else:
                            # Fallback: try original domain or default to Eukarya color for unknown
                            color = DOMAIN_COLORS.get(domain, DOMAIN_COLORS['Eukarya'])
            id_colors.append((id_, color))
            used_colors.add(color)
            if debug:
                print(f"{id_}\tdomain='{domain}'\tphylum='{phylum}'\tspecies={species_clean}\tcolor={color}")
        else:
            # Missing phylum info - default to gray and warn
            if debug:
                print(f"WARNING: {id_} not found in phylum table")
            id_colors.append((id_, "rgba(128,128,128,0.3)"))
    
    # Generate legend based on used colors
    legend_block = generate_legend_block(used_colors)
    
    # Write output file
    with open(color_file, "w") as out_color:
        out_color.write("DATASET_RANGE\n")
        out_color.write("SEPARATOR SPACE\n")
        out_color.write(legend_block)
        out_color.write("\n")
        out_color.write("DATASET_LABEL Example\n")
        out_color.write("DATA\n")
        for id_, color in id_colors:
            out_color.write(f"{id_} {id_} {color}\n")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python dataset-ranges.py input.fasta phylum_table color_file [protein] [mrca_ids_file]")
        sys.exit(1)
    fasta_file = sys.argv[1]
    phylum_table = sys.argv[2]
    color_file = sys.argv[3]
    protein = sys.argv[4] if len(sys.argv) > 4 else None
    mrca_ids_file = sys.argv[5] if len(sys.argv) > 5 else None
    main(fasta_file, phylum_table, color_file, protein, mrca_ids_file=mrca_ids_file)
