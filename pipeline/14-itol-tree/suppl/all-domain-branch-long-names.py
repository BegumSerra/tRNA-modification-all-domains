"""
Generate an iTOL LABELS file with long species names for all leaves on the tree.

For EP-prefixed IDs (eukaryotes), the label shows a shortened ID:
    EP00682_Trypanosoma_cruzi_P003262,EP00682_P003262  Trypanosoma cruzi

For other IDs (WP_, GCA, etc.):
    WP_013558966.1,WP_013558966.1  Anaerolinea thermophila

Reads:
  - iqtree treefile  → extract leaf IDs in tree order
  - phylum table     → map ID → species

Usage:
    python all-domain-branch-long-names.py <treefile> <phylum_table> <output_file>
"""

import re
import sys


def parse_treefile_ids(treefile_path):
    """Extract leaf IDs from a Newick treefile in tree-traversal order."""
    with open(treefile_path) as f:
        tree = f.read().strip()
    leaves = re.findall(r'[\(,]([^():,]+):', tree)
    leaves = [l for l in leaves if not l.replace('.', '').replace('-', '').isdigit()]
    return leaves


def parse_phylum_table(table_path):
    """Parse phylum table: acc → species name."""
    acc_to_species = {}
    with open(table_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                acc = parts[0]
                species = parts[5]
                acc_to_species[acc] = species
    return acc_to_species


def make_short_id(acc):
    """For EP-prefixed IDs, return shortened form: EPxxxxx_Pxxxxxx.
    E.g. EP00682_Trypanosoma_cruzi_P003262 → EP00682_P003262
    For other IDs, return as-is."""
    if acc.startswith("EP"):
        parts = acc.split("_")
        # Last part should be the P-number
        if len(parts) >= 3 and re.match(r'^P\d+$', parts[-1]):
            return f"{parts[0]}_{parts[-1]}"
    return acc


def main(treefile, phylum_table, output_file):
    # 1. Get leaf IDs from treefile (tree order)
    ids = parse_treefile_ids(treefile)

    # 2. Map ID → species from phylum table
    acc_to_species = parse_phylum_table(phylum_table)

    # 3. Write iTOL LABELS file
    matched = 0
    unmatched = 0
    with open(output_file, "w") as out:
        out.write("LABELS\n")
        out.write("SEPARATOR COMMA\n")
        out.write("DATA\n")
        for acc in ids:
            species = acc_to_species.get(acc, "")
            short_id = make_short_id(acc)
            if species:
                # Convert underscores to spaces for display
                species_display = species.replace("_", " ")
                out.write(f"{acc},{short_id}  {species_display}\n")
                matched += 1
            else:
                # No species found — label with just the short ID
                out.write(f"{acc},{short_id}\n")
                unmatched += 1

    print(f"Total leaves: {len(ids)}")
    print(f"Matched to species: {matched}")
    print(f"Unmatched (ID only): {unmatched}")
    print(f"Output saved to: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python all-domain-branch-long-names.py <treefile> <phylum_table> <output_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
