# This script processes the CSV file as described and outputs a space-separated file.
import csv
import sys

# Define allowed localizations and their RGB colors
loc_colors = {
    "Cytoplasm": "rgba(255,0,0,0.5)",         # Red
    "Mitochondrion": "rgba(255,255,0,0.5)",          # Yellow
    "Plastid": "rgba(0,255,0,0.5)",     # Green
    "Nucleus": "rgba(128,128,128,0.5)"           # Gray
}

# Set static thresholds for each localization
thresholds = {
    "Cytoplasm": 0.4695,
    "Nucleus": 0.5275,
    "Mitochondrion": 0.6373,
    "Plastid": 0.6586
}

legend_block = """
LEGEND_TITLE Localization
LEGEND_SCALE 1
LEGEND_POSITION_X 100
LEGEND_POSITION_Y 100
LEGEND_HORIZONTAL 1
LEGEND_SHAPES 1 1 1 1
LEGEND_COLORS rgba(255,0,0,0.5) rgba(255,255,0,0.5) rgba(0,255,0,0.5) rgba(128,128,128,0.5)
LEGEND_LABELS Cytoplasm Mitochondrion Plastid Other
LEGEND_SHAPE_SCALES 1 1 1 1
"""

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, newline='') as csvfile, open(output_file, "w") as outfile:
    outfile.write("DATASET_RANGE\n")
    outfile.write("SEPARATOR SPACE\n")
    outfile.write(legend_block)
    outfile.write("\n")
    outfile.write("DATASET_LABEL Example\n")
    outfile.write("DATA\n")
    reader = csv.DictReader(csvfile)
    cytoplasm_count = 0
    mitochondria_count = 0
    chloroplast_count = 0
    for row in reader:
        protein_id = row["Protein_ID"] if "Protein_ID" in row else row[list(row.keys())[0]]
        
        # Keep only eukaryotic sequences (IDs starting with "EP")
        if not protein_id.startswith("EP"):
            continue
            
        pred_locs = row.get("Localizations", "").replace(",", "|").replace(";", "|").split("|") if "Localizations" in row else []
        pred_locs = [loc.strip() for loc in pred_locs if loc.strip()]
        best_loc = None
        best_diff = float('-inf')
        # Only consider main localizations for threshold comparison
        for loc in pred_locs:
            if loc not in thresholds:
                continue
            score_str = row.get(loc, "")
            if not score_str.strip():
                continue
            try:
                score = float(score_str)
                threshold = thresholds[loc]
                diff = score - threshold
                if diff > best_diff:
                    best_diff = diff
                    best_loc = loc
            except (ValueError, KeyError):
                continue
        if best_loc:
            color = loc_colors.get(best_loc, "rgba(128,128,128,0.5)")
            outfile.write(f"{protein_id} {protein_id} {color}\n")
            if best_loc.lower() == "cytoplasm":
                cytoplasm_count += 1
            elif best_loc.lower() == "mitochondrion":
                mitochondria_count += 1
            elif best_loc.lower() == "plastid":
                chloroplast_count += 1
        else:
            # If none of the main localizations are above threshold, assign as Other (gray)
            outfile.write(f"{protein_id} {protein_id} rgba(128,128,128,0.5)\n")
    # Write summary table to a separate file
    summary_file = outfile.name + ".summary.tsv"
    with open(summary_file, "w") as sumf:
        sumf.write("Cytoplasm\tMitochondria\tChloroplast\n")
        sumf.write(f"{cytoplasm_count}\t{mitochondria_count}\t{chloroplast_count}\n")
