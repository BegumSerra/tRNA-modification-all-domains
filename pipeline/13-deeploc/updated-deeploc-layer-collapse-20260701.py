#!/usr/bin/env python3
import sys
import csv

# Color definitions for each group
GROUPS = {
    "Cytoplasm/Nucleus":    "#E2D0F1",    #(combines Cytoplasm, Nucleus, and Cytoplasm|Nucleus)
    "Mitochondrion":        "#ff0097",      # red-orange
    "Plastid":              "#00aba9",     # lime green
    "Other":                "#D4F5DD",   # C2B7AB,   
}


def classify(localizations_str):
    """Classify a localization string into one of the defined groups"""
    locs = set(l.strip() for l in localizations_str.split("|") if l.strip())

    has_cyto  = "Cytoplasm"     in locs
    has_nuc   = "Nucleus"       in locs
    has_mito  = "Mitochondrion" in locs
    has_plast = "Plastid"       in locs

    # Group Cytoplasm, Nucleus, and Cytoplasm|Nucleus together
    if has_cyto or has_nuc:
        return "Cytoplasm/Nucleus"
    elif has_mito:
        return "Mitochondrion"
    elif has_plast:
        return "Plastid"
    else:
        return "Other"


def main(deeploc_csv, output_file):
    protein_groups = {}
    other_details  = {}

    with open(deeploc_csv, newline="") as f:
        for row in csv.DictReader(f):
            protein_id = row.get("Protein_ID", "").strip()
            locs       = row.get("Localizations", "").strip()
            if not protein_id:
                continue
            group = classify(locs)
            protein_groups[protein_id] = group
            if group == "Other":
                other_details[protein_id] = locs

    # Stats
    from collections import Counter
    counts = Counter(protein_groups.values())
    print(f"Total proteins:     {len(protein_groups)}")
    for g, c in counts.most_common():
        print(f"  {g:<25} {c}")

    if other_details:
        print(f"\nProteins in 'Other' group:")
        for pid, locs in other_details.items():
            print(f"  {pid} → {locs}")

    # Legend — only groups that appear in data
    legend_groups = [g for g in GROUPS if any(v == g for v in protein_groups.values())]

    shapes  = " ".join(["1"] * len(legend_groups))
    colors  = " ".join([GROUPS[g] for g in legend_groups])
    labels  = " ".join([g.replace("|", "_").replace(" ", "_") for g in legend_groups])

    with open(output_file, "w") as out:
        out.write("DATASET_RANGE\n")
        out.write("SEPARATOR SPACE\n")
        out.write("DATASET_LABEL Localization\n\n")
        out.write("LEGEND_TITLE Localization\n")
        out.write("LEGEND_SCALE 1\n")
        out.write("LEGEND_POSITION_X 100\n")
        out.write("LEGEND_POSITION_Y 100\n")
        out.write("LEGEND_HORIZONTAL 1\n")
        out.write(f"LEGEND_SHAPES {shapes}\n")
        out.write(f"LEGEND_COLORS {colors}\n")
        out.write(f"LEGEND_LABELS {labels}\n")
        out.write(f"LEGEND_SHAPE_SCALES {' '.join(['1'] * len(legend_groups))}\n\n")
        out.write("DATA\n")
        out.write("#ID LABEL COLOR\n")

        for protein_id, group in protein_groups.items():
            color = GROUPS[group]
            out.write(f"{protein_id} {protein_id} {color}\n")

    print(f"\nOutput written to: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python deeploc_itol.py <deeploc_csv> <output>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])