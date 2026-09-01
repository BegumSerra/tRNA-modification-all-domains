# This script processes the CSV file as described and outputs a space-separated file.
import csv
import sys

# Define allowed localizations and their RGB colors
loc_colors = {
    "Cytoplasm": "rgba(255,0,0,0.5)",         # Red
    "Mitochondrion": "rgba(255,255,0,0.5)",          # Yellow
    "Plastid": "rgba(0,255,0,0.5)",     # Green
    "Nucleus": "rgba(0,0,255,0.5)"           # Blue
}

legend_block = """
LEGEND_TITLE Localization
LEGEND_SCALE 1
LEGEND_POSITION_X 100
LEGEND_POSITION_Y 100
LEGEND_HORIZONTAL 1
LEGEND_SHAPES 1 1 1 1
LEGEND_COLORS rgba(255,0,0,0.5) rgba(255,255,0,0.5) rgba(0,255,0,0.5) rgba(0,0,255,0.5)
LEGEND_LABELS Cytoplasm Mitochondrion Plastid Nucleus
LEGEND_SHAPE_SCALES 1 1 1 1
"""

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, newline='') as csvfile, open(output_file, "w") as outfile:
    # Write metadata lines at the top of the output file
    outfile.write("DATASET_RANGE\n")
    outfile.write("SEPARATOR SPACE\n")
    outfile.write(legend_block)
    outfile.write("\n")
    outfile.write("DATASET_LABEL Example\n")
    outfile.write("DATA\n")
    reader = csv.reader(csvfile)
    header = next(reader)  # skip header
    for row in reader:
        if len(row) < 2:
            continue
        protein_id = row[0]
        localizations = row[1].split("|")
        filtered = [loc for loc in localizations if loc in loc_colors]
        # Only keep rows where there is exactly one allowed localization
        if len(filtered) == 1:
            loc = filtered[0]
            color = loc_colors[loc]
            outfile.write(f"{protein_id} {protein_id} {color}\n")

