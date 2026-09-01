"""
Generate an iTOL DATASET_TEXT file that places a text label for EVERY
consecutive run of a group — not just the largest one.

If a group appears in multiple separate places on the tree, each
occurrence gets its own label at the middle of that run.

Works for any single domain — the domain-specific gradient file already
contains only the IDs for that domain.

Usage:
    python all-domain-groups-labels.py <gradient_file> <domain_label> <output_file>

    domain_label:  e.g. "Archaea", "Eukaryota", "Bacteria"
                   (used only in legend titles / dataset labels)
"""

import sys


def parse_gradient_file(gradient_path):
    """Parse the gradient range file to get ID → gradient_color mapping and legend info."""
    id_color = {}
    legend_colors = []
    legend_labels = []

    in_data = False
    with open(gradient_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("LEGEND_COLORS"):
                legend_colors = line.split(None, 1)[1].split()
            elif line.startswith("LEGEND_LABELS"):
                legend_labels = line.split(None, 1)[1].split()
            elif line == "DATA":
                in_data = True
                continue
            elif in_data:
                # Format: ID LABEL COLOR GRADIENT_COLOR
                parts = line.split()
                if len(parts) >= 4:
                    acc = parts[0]
                    gradient_color = parts[3]
                    id_color[acc] = gradient_color

    return id_color, legend_colors, legend_labels


def main(gradient_file, domain_label, output_file):
    # Parse gradient file (already domain-specific)
    id_color, legend_colors, legend_labels = parse_gradient_file(gradient_file)

    if not id_color:
        # Empty gradient file — write placeholder
        with open(output_file, "w") as out:
            out.write(f"# No {domain_label} groups to label\n")
        print(f"No IDs in gradient file — empty placeholder written → {output_file}")
        return

    # Build color → group label from legend
    color_to_group = {}
    for c, l in zip(legend_colors, legend_labels):
        color_to_group[c] = l.replace("_", " ")

    # Map each ID to its group name via its gradient color
    id_to_group = {}
    for acc, color in id_color.items():
        group = color_to_group.get(color)
        if group:
            id_to_group[acc] = group

    print(f"Total IDs in gradient file: {len(id_color)}")
    print(f"Groups in legend: {len(legend_labels)}")

    # Write DATASET_TEXT iTOL file
    with open(output_file, "w") as out:
        out.write("DATASET_TEXT\n")
        out.write(f"#Text labels for {domain_label} groups, colored to match the gradient ranges\n")
        out.write("SEPARATOR SPACE\n\n")

        # Legend
        out.write(f"LEGEND_TITLE {domain_label}_Groups\n")
        out.write("LEGEND_SHAPES " + " ".join(["1"] * len(legend_colors)) + "\n")
        out.write("LEGEND_COLORS " + " ".join(legend_colors) + "\n")
        out.write("LEGEND_LABELS " + " ".join(legend_labels) + "\n")
        out.write("LEGEND_SHAPE_SCALES " + " ".join(["1"] * len(legend_colors)) + "\n\n")

        # Dataset settings
        out.write(f"DATASET_LABEL {domain_label}_Group_Labels\n")
        out.write("COLOR #000000\n\n")

        # Display settings
        out.write("MARGIN 5\n")
        out.write("SHOW_INTERNAL 0\n")
        out.write("ALIGN_TO_TREE 1\n")
        out.write("SIZE_FACTOR 1\n")
        out.write("STRAIGHT_LABELS 1\n")
        out.write("LABELS_BELOW 0\n\n")

        # Data section
        out.write("DATA\n")
        out.write("#ID label position color style size_factor rotation\n")

        # Build ALL consecutive runs (tree order from gradient file)
        gradient_ids = list(id_color.keys())  # preserves insertion order = tree order

        runs = []
        for acc in gradient_ids:
            if acc in id_to_group:
                group_name = id_to_group[acc]
                if runs and runs[-1][0] == group_name:
                    runs[-1][1].append(acc)
                else:
                    runs.append((group_name, [acc]))

        # Write a label for EVERY consecutive run (even if a group appears multiple times)
        labels_written = 0
        from collections import Counter
        run_counter = Counter()
        for group_name, accs in runs:
            run_counter[group_name] += 1
            mid = len(accs) // 2
            rep_acc = accs[mid]
            color = id_color[rep_acc]
            occurrence = run_counter[group_name]
            print(f"  {group_name} (run {occurrence}): {len(accs)} IDs, representative={rep_acc}")
            out.write(f"{rep_acc} {group_name} -1 {color} bold 1 0\n")
            labels_written += 1

    # Count groups that appear more than once
    multi = {g: c for g, c in run_counter.items() if c > 1}
    unique_groups = len(run_counter)

    print(f"\nOutput saved to: {output_file}")
    print(f"Unique {domain_label} groups: {unique_groups}")
    print(f"Total labels written (all runs): {labels_written}")
    if multi:
        print(f"Groups with multiple runs:")
        for g, c in sorted(multi.items()):
            print(f"  {g}: {c} runs")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python all-domain-groups-labels.py <gradient_file> <domain_label> <output_file>")
        print("  gradient_file: domain-specific gradient iTOL file")
        print("  domain_label:  Archaea, Eukaryota, or Bacteria")
        print("  output_file:   output iTOL DATASET_TEXT file")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
