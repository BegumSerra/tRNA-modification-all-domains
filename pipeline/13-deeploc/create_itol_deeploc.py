#!/usr/bin/env python3
"""
Convert DeepLoc CSV results to iTOL dataset_range format
Usage: python create_itol_deeploc.py input.csv output.txt
"""

import csv
import sys
import argparse

# Define localization colors - 6 groups
LOCALIZATION_COLORS = {
    "Cytoplasm": "rgba(255,0,0,0.5)",            # Red - Cytoplasm only
    "Nucleus": "rgba(0,0,255,0.5)",              # Blue - Nucleus only
    "Cytoplasm|Nucleus": "rgba(128,0,128,0.5)",  # Purple - Both Cytoplasm and Nucleus
    "Mitochondrion": "rgba(255,255,0,0.5)",      # Yellow - Mitochondrion
    "Plastid": "rgba(0,255,0,0.5)",              # Green - Plastid
    "Other": "rgba(128,128,128,0.5)"             # Gray - Everything else
}


def parse_localization(row):
    """
    Determine the localization for a protein based on Localizations column.
    Maps to 6 distinct groups:
    1. Cytoplasm only
    2. Nucleus only
    3. Cytoplasm|Nucleus (dual prediction)
    4. Mitochondrion
    5. Plastid
    6. Other
    
    Args:
        row: CSV row dictionary
    
    Returns:
        tuple: (localization_name, color)
    """
    # Get predicted localizations from "Localizations" column
    localizations = row.get("Localizations", "").strip()
    
    if not localizations:
        return "Other", LOCALIZATION_COLORS["Other"]
    
    # Normalize separators to |
    localizations = localizations.replace(",", "|").replace(";", "|")
    pred_locs = [loc.strip() for loc in localizations.split("|") if loc.strip()]
    
    # Check for exact matches first
    # Cytoplasm|Nucleus (dual prediction)
    if set(pred_locs) == {"Cytoplasm", "Nucleus"}:
        return "Cytoplasm|Nucleus", LOCALIZATION_COLORS["Cytoplasm|Nucleus"]
    
    # Single Cytoplasm
    if pred_locs == ["Cytoplasm"]:
        return "Cytoplasm", LOCALIZATION_COLORS["Cytoplasm"]
    
    # Single Nucleus
    if pred_locs == ["Nucleus"]:
        return "Nucleus", LOCALIZATION_COLORS["Nucleus"]
    
    # Mitochondrion (with or without other predictions - prioritize organelle)
    if "Mitochondrion" in pred_locs:
        return "Mitochondrion", LOCALIZATION_COLORS["Mitochondrion"]
    
    # Plastid (with or without other predictions - prioritize organelle)
    if "Plastid" in pred_locs:
        return "Plastid", LOCALIZATION_COLORS["Plastid"]
    
    # Everything else goes to Other
    return "Other", LOCALIZATION_COLORS["Other"]


def create_itol_file(input_csv, output_txt, protein_name="Protein", 
                     eukaryotes_only=True):
    """
    Convert DeepLoc CSV to iTOL dataset_range format.
    
    Args:
        input_csv: Path to DeepLoc CSV file
        output_txt: Path to output iTOL file
        protein_name: Name for the dataset label
        eukaryotes_only: If True, only include proteins starting with "EP"
    """
    
    # Count localizations
    loc_counts = {loc: 0 for loc in LOCALIZATION_COLORS.keys()}
    
    with open(input_csv, 'r', newline='') as csvfile, open(output_txt, 'w') as outfile:
        # Write iTOL header
        outfile.write("DATASET_RANGE\n")
        outfile.write("SEPARATOR SPACE\n\n")
        
        # Write legend
        outfile.write("LEGEND_TITLE Localization\n")
        outfile.write("LEGEND_SCALE 1\n")
        outfile.write("LEGEND_POSITION_X 100\n")
        outfile.write("LEGEND_POSITION_Y 100\n")
        outfile.write("LEGEND_HORIZONTAL 1\n")
        
        # Legend with 6 groups
        legend_labels = ["Cytoplasm", "Nucleus", "Cytoplasm|Nucleus", "Mitochondrion", "Plastid", "Other"]
        legend_colors = [
            "rgba(255,0,0,0.5)",       # Red
            "rgba(0,0,255,0.5)",       # Blue
            "rgba(128,0,128,0.5)",     # Purple
            "rgba(255,255,0,0.5)",     # Yellow
            "rgba(0,255,0,0.5)",       # Green
            "rgba(128,128,128,0.5)"    # Gray
        ]
        
        outfile.write(f"LEGEND_SHAPES {' '.join(['1'] * 6)}\n")
        outfile.write(f"LEGEND_COLORS {' '.join(legend_colors)}\n")
        outfile.write(f"LEGEND_LABELS {' '.join(legend_labels)}\n")
        outfile.write(f"LEGEND_SHAPE_SCALES {' '.join(['1'] * 6)}\n\n")
        
        outfile.write(f"DATASET_LABEL {protein_name}\n")
        outfile.write("DATA\n")
        
        # Process CSV rows
        reader = csv.DictReader(csvfile)
        for row in reader:
            protein_id = row.get("Protein_ID", "")
            
            # Filter to eukaryotes if requested
            if eukaryotes_only and not protein_id.startswith("EP"):
                continue
            
            # Get localization and color
            loc, color = parse_localization(row)
            loc_counts[loc] += 1
            
            # Write iTOL line: protein_id protein_id color
            outfile.write(f"{protein_id} {protein_id} {color}\n")
    
    # Print summary
    print(f"Created iTOL file: {output_txt}")
    print(f"Total proteins: {sum(loc_counts.values())}")
    print("\nLocalization counts:")
    for loc, count in sorted(loc_counts.items(), key=lambda x: -x[1]):
        if count > 0:
            print(f"  {loc}: {count}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert DeepLoc CSV to iTOL dataset_range format"
    )
    parser.add_argument("input_csv", help="Input DeepLoc CSV file")
    parser.add_argument("output_txt", help="Output iTOL file")
    parser.add_argument("--protein-name", default="Protein", 
                       help="Dataset label name (default: Protein)")
    parser.add_argument("--all-proteins", action="store_true",
                       help="Include all proteins (not just eukaryotes starting with EP)")
    
    args = parser.parse_args()
    
    create_itol_file(
        args.input_csv,
        args.output_txt,
        protein_name=args.protein_name,
        eukaryotes_only=not args.all_proteins
    )


if __name__ == "__main__":
    main()
