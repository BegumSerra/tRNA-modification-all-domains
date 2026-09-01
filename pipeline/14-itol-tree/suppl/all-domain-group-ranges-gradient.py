"""
Generate per-domain iTOL DATASET_RANGE (gradient) files.

Produces up to three separate files — one each for Archaea, Eukaryota,
and Bacteria — but only if that domain has sequences on the tree.

Reads:
  - iqtree treefile  → extract leaf IDs in tree order
  - phylum table     → map ID → (domain, species)
  - color.tsv        → map species → (Group, Colorcode)

Usage:
    python archaea-order-ranges-gradient.py <treefile> <phylum_table> <color_tsv> \\
           <archaea_output> <eukaryota_output> <bacteria_output>
"""

import re
import sys


# Map phylum-table domain names to short tags used in filenames / messages
DOMAIN_MAP = {
    "Archaea": "archaea",
    "Eukaryota": "eukaryota",
    "Bacteria": "bacteria",
}


def parse_treefile_ids(treefile_path):
    """Extract leaf IDs from a Newick treefile in tree-traversal order."""
    with open(treefile_path) as f:
        tree = f.read().strip()
    leaves = re.findall(r'[\(,]([^():,]+):', tree)
    leaves = [l for l in leaves if not l.replace('.', '').replace('-', '').isdigit()]
    return leaves


def parse_phylum_table(table_path):
    """Parse phylum table: acc → (domain, species)."""
    acc_info = {}
    with open(table_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                acc = parts[0]
                domain = parts[1]
                species = parts[5]
                acc_info[acc] = (domain, species)
    return acc_info


def parse_color_table(color_path):
    """Parse color.tsv: species → (Group, Colorcode).

    Returns:
        species_to_group: dict  species_name → (group, colorcode)
        group_colors: dict      group_name → colorcode  (first occurrence wins)
        group_order: list       groups in file order (unique, preserves first appearance)
    """
    species_to_group = {}
    group_colors = {}
    group_order = []
    with open(color_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                species = parts[4]
                group = parts[5]
                colorcode = parts[6]
                species_to_group[species] = (group, colorcode)
                if group not in group_colors:
                    group_colors[group] = colorcode
                    group_order.append(group)
    return species_to_group, group_colors, group_order


def generate_legend_block(group_colors, group_order, domain_label):
    """Generate iTOL legend block from group→color mapping."""
    groups = [g for g in group_order if g in group_colors]
    colors = [group_colors[g] for g in groups]
    shapes = " ".join(["1"] * len(groups))
    colors_str = " ".join(colors)
    labels_str = " ".join([g.replace(" ", "_") for g in groups])
    scales = " ".join(["1"] * len(groups))

    legend = f"""DATASET_RANGE
SEPARATOR SPACE
LEGEND_TITLE {domain_label}_Groups
LEGEND_SCALE 1
LEGEND_POSITION_X 100
LEGEND_POSITION_Y 100
LEGEND_HORIZONTAL 1
LEGEND_SHAPES {shapes}
LEGEND_COLORS {colors_str}
LEGEND_LABELS {labels_str}
LEGEND_SHAPE_SCALES {scales}

DATASET_LABEL {domain_label}_Groups
DATA
#ID LABEL COLOR GRADIENT_COLOR
"""
    return legend


def write_domain_gradient(ids, id_to_group, id_to_color, domain_ids,
                          group_colors, group_order, domain_label, output_file):
    """Write a single-domain gradient file. Returns number of IDs written."""
    domain_set = set(domain_ids)

    # Collect groups present for this domain (in color.tsv order)
    seen = set()
    for acc in ids:
        if acc in domain_set and acc in id_to_group:
            seen.add(id_to_group[acc])
    legend_groups = [g for g in group_order if g in seen]

    with open(output_file, "w") as out:
        out.write(generate_legend_block(group_colors, legend_groups, domain_label))
        for acc in ids:
            if acc in domain_set and acc in id_to_color:
                color = id_to_color[acc]
                out.write(f"{acc} {acc} #ffffff {color}\n")

    return len([a for a in ids if a in domain_set and a in id_to_color]), len(legend_groups)


def main(treefile, phylum_table, color_tsv, arc_out, euk_out, bact_out):
    # 1. Get leaf IDs from treefile (tree order)
    ids = parse_treefile_ids(treefile)

    # 2. Map ID → (domain, species) from phylum table
    acc_info = parse_phylum_table(phylum_table)

    # 3. Map species → (Group, Colorcode) from color.tsv
    species_to_group, group_colors, group_order = parse_color_table(color_tsv)

    # 4. Map each leaf ID → group + color + domain
    id_to_group = {}
    id_to_color = {}
    id_to_domain = {}
    unmatched = []
    for acc in ids:
        info = acc_info.get(acc)
        if info is None:
            unmatched.append((acc, "not in phylum table"))
            continue
        domain, species = info
        id_to_domain[acc] = domain
        match = species_to_group.get(species) or species_to_group.get(species.replace("_", " "))
        if match:
            group, color = match
            id_to_group[acc] = group
            id_to_color[acc] = color
        else:
            unmatched.append((acc, species))

    # 5. Split IDs by domain
    domain_ids = {"Archaea": [], "Eukaryota": [], "Bacteria": []}
    for acc in ids:
        d = id_to_domain.get(acc)
        if d in domain_ids:
            domain_ids[d].append(acc)

    print(f"Total leaves in tree: {len(ids)}")
    print(f"Matched to group: {len(id_to_group)}")
    print(f"Unmatched: {len(unmatched)}")
    if unmatched[:5]:
        print("Sample unmatched:")
        for acc, sp in unmatched[:5]:
            print(f"  {acc}: {sp}")

    # 6. Write per-domain gradient files (only if domain has sequences)
    outputs = [
        ("Archaea", "Archaea", arc_out),
        ("Eukaryota", "Eukaryota", euk_out),
        ("Bacteria", "Bacteria", bact_out),
    ]

    for domain_key, domain_label, out_path in outputs:
        d_ids = domain_ids[domain_key]
        if d_ids:
            n_ids, n_groups = write_domain_gradient(
                ids, id_to_group, id_to_color, d_ids,
                group_colors, group_order, domain_label, out_path
            )
            print(f"\n{domain_label}: {len(d_ids)} leaves, {n_ids} colored, {n_groups} groups → {out_path}")
        else:
            # Write an empty placeholder so snakemake doesn't fail
            with open(out_path, "w") as out:
                out.write(f"# No {domain_label} sequences on this tree\n")
            print(f"\n{domain_label}: 0 leaves — empty placeholder written → {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python archaea-order-ranges-gradient.py <treefile> <phylum_table> <color_tsv> "
              "<archaea_out> <eukaryota_out> <bacteria_out>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
