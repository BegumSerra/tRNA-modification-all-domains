import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
import sys

# Input file
input_file = sys.argv[1]
protein_group_file = sys.argv[2]
output_file = sys.argv[3]

# Load the CSV file
df_large = pd.read_csv(input_file, header=0, sep="\t")

# Define domain colors
domain_colors = {
    'Archaea': '#008CE3',   # Blue
    'Eukarya': '#FF6A38',   # Orange
    'Bacteria': '#CE0055'   # Pink
}

# Light colors for half-filled circles (50% lighter)
domain_colors_light = {
    'Archaea': '#80C6F1',   # Light blue
    'Eukarya': '#FFB59C',   # Light orange
    'Bacteria': '#E680AA'   # Light pink
}

# Create a mapping from group to domain (use first occurrence)
group_to_domain = df_large.groupby('Group')['Domain'].first().to_dict()

# Ensure the required columns exist
if 'Group' not in df_large.columns:
    raise ValueError("The CSV file must have a 'Group' column.")
if 'Domain' not in df_large.columns:
    raise ValueError("The CSV file must have a 'Domain' column.")

# Remove rows with missing Group values
df_large = df_large.dropna(subset=['Group'])
df_large = df_large[df_large['Group'] != 'Unknown']

# Extract protein names (all columns except taxonomy columns)
taxonomy_cols = ['Domain', 'Group', 'Phylum', 'Order', 'Family', 'Species']
proteins_to_remove = []
protein_names = [col for col in df_large.columns if col not in taxonomy_cols and col not in proteins_to_remove]

# Load protein functional groups from file
# Each line is a comma-separated group of proteins
protein_func_groups = []
with open(protein_group_file, 'r') as f:
    for line in f:
        proteins = [p.strip() for p in line.strip().rstrip(',').split(',') if p.strip()]
        if proteins:
            protein_func_groups.append(proteins)

# Build ordered protein list from the group file (only proteins that exist in the data and have non-zero presence)
all_proteins_in_groups = []
for func_group in protein_func_groups:
    for p in func_group:
        all_proteins_in_groups.append(p)

# Filter to valid proteins (exist in data and have non-zero presence)
valid_proteins = []
for p in all_proteins_in_groups:
    if p in protein_names and df_large[p].sum() > 0:
        valid_proteins.append(p)

# Add remaining proteins not in the group file
for p in protein_names:
    if p not in valid_proteins and df_large[p].sum() > 0:
        valid_proteins.append(p)

# Build a mapping: protein -> functional group index
protein_to_func_group = {}
for grp_idx, func_group in enumerate(protein_func_groups):
    for p in func_group:
        protein_to_func_group[p] = grp_idx

print(f"Protein functional groups loaded: {len(protein_func_groups)} groups")
print(f"Valid proteins: {len(valid_proteins)}")
print(f"Protein order: {valid_proteins}")

# Step 2: Filter the dataframe to only include the valid proteins
df_filtered = df_large[taxonomy_cols + valid_proteins]

# Step 3: Count species per group
species_per_group = df_filtered.groupby("Group").size()

# Step 4: Group by Group and calculate the proportion of presence for each protein
df_grouped_large = df_filtered.groupby("Group")[valid_proteins].apply(lambda x: (x > 0).sum() / len(x))

# Explicitly reorder columns to match the ordered protein list
df_grouped_large = df_grouped_large[valid_proteins]

print(f"Protein order in plot: {list(df_grouped_large.columns)}")

# Step 5: Preserve the group order as it appears in the data
unique_group_order = list(df_filtered['Group'].unique())
unique_group_order = [g for g in unique_group_order if g in df_grouped_large.index]
df_grouped_large = df_grouped_large.loc[unique_group_order]

# Step 6: Create group labels with species count in parentheses
group_labels_with_count = [f"{group} ({species_per_group[group]})" for group in df_grouped_large.index]

# Step 7: Create the figure and axis
fig, ax = plt.subplots(figsize=(15, 10))

# Define circle and arc sizes
circle_radius = 0.3
domain_gap = 1.5  # Extra space between different domains
protein_group_gap = 0.8  # Extra space between protein functional groups

# Step 8: Calculate y-positions with gaps between domains
y_positions = []
y = 0
prev_domain = None
for group in df_grouped_large.index:
    domain = group_to_domain.get(group, 'Bacteria')
    if prev_domain is not None and domain != prev_domain:
        y += domain_gap  # Add extra space between domains
    y_positions.append(y)
    prev_domain = domain
    y += 1

# Step 8b: Calculate x-positions with gaps between protein functional groups
x_positions = []
x = 0
prev_func_grp = protein_to_func_group.get(df_grouped_large.columns[0], -1)
for j, protein in enumerate(df_grouped_large.columns):
    func_grp = protein_to_func_group.get(protein, -1)
    if j > 0 and func_grp != prev_func_grp:
        x += protein_group_gap  # Add gap between functional groups
    x_positions.append(x)
    prev_func_grp = func_grp
    x += 1

print(f"X-positions: {x_positions[:10]}... (total {len(x_positions)})")

# Step 9: Iterate over groups (rows) and proteins (columns)
for i, group in enumerate(df_grouped_large.index):
    # Get domain for this group
    domain = group_to_domain.get(group, 'Bacteria')
    fill_color = domain_colors.get(domain, '#CE0055')
    light_color = domain_colors_light.get(domain, '#E680AA')
    
    for j, protein in enumerate(df_grouped_large.columns):
        proportion = df_grouped_large.loc[group, protein]
        x = x_positions[j]  # Protein position (column) with functional group gaps
        y = y_positions[i]  # Group position (row) with domain gaps

        # Full circle for ≥ 50% presence
        if proportion >= 0.5:
            circle = Circle((x, y), circle_radius, edgecolor='lightgray', facecolor=fill_color, alpha=0.8)
            ax.add_patch(circle)
        # Semi-circle for < 50% presence
        elif proportion > 0:
            half_filled_circle = Wedge((x, y), circle_radius, 0, 180, color=light_color, alpha=0.8)
            ax.add_patch(half_filled_circle)
            outline_circle = Circle((x, y), circle_radius, edgecolor='lightgray', facecolor='none', linewidth=1)
            ax.add_patch(outline_circle)
        # Outline circle for absent proteins (0 presence)
        else:
            outline_circle = Circle((x, y), circle_radius, edgecolor='lightgray', facecolor='none', linewidth=1)
            ax.add_patch(outline_circle)

# Step 10: Set labels with species count in parentheses
ax.set_xticks(x_positions)
ax.set_xticklabels(df_grouped_large.columns, rotation=90)
ax.set_yticks(y_positions)
ax.set_yticklabels(group_labels_with_count)

# Color the y-axis labels based on domain
for i, (tick, group) in enumerate(zip(ax.get_yticklabels(), df_grouped_large.index)):
    domain = group_to_domain.get(group, 'Bacteria')
    tick.set_color(domain_colors.get(domain, '#CE0055'))

# Step 11: Adjust axis limits and aspect ratio
ax.set_xlim(-1, max(x_positions) + 1)
ax.set_ylim(-1, max(y_positions) + 1)
ax.invert_yaxis()  # Invert y-axis so first group in table appears at the top
ax.set_aspect('equal')

# Print the number of proteins present in each group
print("Number of proteins present in each group:")
print(df_filtered.groupby('Group')[valid_proteins].apply(lambda x: (x > 0).sum(axis=1).sum()))

# Show the plot
plt.tight_layout()
# Save the figure
plt.savefig(output_file, format="pdf", dpi=300)
print(f"Plot saved to {output_file}")
plt.close()