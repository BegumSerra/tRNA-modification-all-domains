from Bio import Phylo
import sys
import re
import os

TREEFILE = sys.argv[1]
output_file = sys.argv[2]
protein = sys.argv[3] if len(sys.argv) > 3 else None

# Automatically detect protein from input tree filename if not provided
if protein is None:
    filename = os.path.basename(TREEFILE).lower()
    if 'mnmg' in filename:
        protein = 'MnmG'
    elif 'mnme' in filename:
        protein = 'MnmE'
    elif 'tils' in filename:
        protein = 'TilS'
    else:
        protein = None
    print(f"[Auto-detected protein argument: {protein}]")

# Define target tips for each protein
MNMG_TIPS = [
    "WP_164928623.1",
    "WP_099799201.1",
    "WP_041566549.1",
    "WP_254035056.1",
    "WP_011821156.1"
]
MNME_TIPS = [
    "WP_011142447.1",
    "WP_099799239.1",
    "WP_012410642.1",
    "WP_227350184.1",
    "WP_011819541.1"
]

TILS_TIPS = [
    "WP_011821133.1",
    "WP_099798688.1",
    "WP_164929556.1",
    "WP_012412192.1",
    "WP_254035740.1"
]

def get_target_tips(protein):
    if protein == "MnmG":
        return MNMG_TIPS
    elif protein == "MnmE":
        return MNME_TIPS
    elif protein == "TilS":
        return TILS_TIPS
    else:
        return []

TARGET_TIPS = get_target_tips(protein)
print(f"Target tips for this protein: {TARGET_TIPS}")

# Read tree
tree = Phylo.read(TREEFILE, "newick")

# Print all tip names in the tree for debugging
print("All tip names in the tree:")
for term in tree.get_terminals():
    print(term.name)

# Write all tip names to a file for inspection
with open("all_tree_tips.txt", "w") as tipfile:
    for term in tree.get_terminals():
        tipfile.write(term.name + "\n")

# Only keep tips that are present in the tree
all_tree_tips = set(term.name for term in tree.get_terminals())
present_tips = [tip for tip in TARGET_TIPS if tip in all_tree_tips]
missing_tips = [tip for tip in TARGET_TIPS if tip not in all_tree_tips]
print("\nTarget tips present in tree:", present_tips)
print("Target tips missing from tree:", missing_tips)

if len(present_tips) < 2:
    print("Not enough target tips present in the tree to find MRCA.")
    with open(output_file, 'w') as f:
        pass  # Write empty file
    sys.exit(0)

# Find MRCA
mrca = tree.common_ancestor(present_tips)
print("MRCA clade:", mrca)
print("MRCA tips:", [t.name for t in mrca.get_terminals()])

def find_parent(tree, node):
    for clade in tree.find_clades(order="level"):
        if node in clade.clades:
            return clade
    return None

parent = find_parent(tree.root, mrca)
if parent is None:
    print("MRCA is the root; no sister clade.")
    with open(output_file, 'w') as f:
        pass
    sys.exit(0)

# Find sister clade (the other child of the parent)
sister = None
for child in parent.clades:
    if child is not mrca:
        sister = child
        break

if sister is None:
    print("No true sister clade found.")
    with open(output_file, 'w') as f:
        pass
    sys.exit(0)
else:
    sister_tips = [term.name for term in sister.get_terminals()]
    print("Sister clade tips:", sister_tips)
    mrca_tips = [t.name for t in mrca.get_terminals()]
    final_list = list(set(mrca_tips) - set(sister_tips))
    print("Final MRCA tips (excluding sister clade):", final_list)

    # Filter final_list to only include tips starting with 'EP' followed by numbers
    filtered_final_list = [tip for tip in final_list if re.match(r'^EP\d+', tip)]
    print("Filtered final MRCA tips (EP only):", filtered_final_list)

    # Write the filtered list to the output file, one ID per line
    with open(output_file, 'w') as f:
        for tip in filtered_final_list:
            f.write(f"{tip}\n")
