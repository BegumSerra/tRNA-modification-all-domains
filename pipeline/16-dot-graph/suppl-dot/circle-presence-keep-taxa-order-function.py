import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
import sys

# Input file
input_file = sys.argv[1]
species_order_file = sys.argv[2]
output_file = sys.argv[3]

# Load the CSV file
df_large = pd.read_csv(input_file, header=0, sep = ";")

# Load phylum order file (one phylum per line)
with open(species_order_file, 'r') as f:
    phylum_order_from_file = [line.strip() for line in f if line.strip()]

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

# Create a mapping from phylum to domain
phylum_to_domain = df_large.groupby('Phylum')['Domain'].first().to_dict()
# Ensure the required columns exist
if 'Phylum' not in df_large.columns:
    raise ValueError("The CSV file must have a 'Phylum' column.")
if 'Domain' not in df_large.columns:
    raise ValueError("The CSV file must have a 'Domain' column.")

# Handle missing values in the 'Phylum' column: replace NaN with 'Unknown'
df_large['Phylum'] = df_large['Phylum'].fillna('Unknown')

# Extract protein names (all columns except 'Domain' and 'Phylum')
protein_names = [col for col in df_large.columns if col not in ['Domain', 'Phylum']]

# Step 1: Create a list of proteins with a non-zero presence across all species
valid_proteins = [protein for protein in protein_names if df_large[protein].sum() > 0]

# Step 2: Filter the dataframe to only include the valid proteins
df_filtered = df_large[['Domain', 'Phylum'] + valid_proteins]

# Update the protein names list after filtering
protein_names_filtered = [col for col in df_filtered.columns if col not in ['Domain', 'Phylum']]

# Step 3: Count species per phylum
species_per_phylum = df_filtered.groupby("Phylum").size()

# Step 4: Group by phylum and calculate the proportion of presence for each protein
# This step does not sort the phylum order but keeps it as in the input file
df_grouped_large = df_filtered.groupby("Phylum")[protein_names_filtered].apply(lambda x: (x > 0).sum() / len(x))

# Step 5: Preserve the phylum order from species_order_eba.txt
# Filter to only include phyla that exist in the data
unique_phyla_order = [p for p in phylum_order_from_file if p in df_grouped_large.index]
df_grouped_large = df_grouped_large.loc[unique_phyla_order]

# Step 6: Create phylum labels with species count in parentheses
phylum_labels_with_count = [f"{phylum} ({species_per_phylum[phylum]})" for phylum in df_grouped_large.index]

# Step 7: Create the figure and axis
fig, ax = plt.subplots(figsize=(15, 10))

# Define circle and arc sizes
circle_radius = 0.3

# Step 6: Iterate over phyla (rows) and proteins (columns)
for i, phylum in enumerate(df_grouped_large.index):
    # Get domain for this phylum
    domain = phylum_to_domain.get(phylum, 'Bacteria')
    fill_color = domain_colors.get(domain, '#CE0055')
    light_color = domain_colors_light.get(domain, '#E680AA')
    
    for j, protein in enumerate(df_grouped_large.columns):
        proportion = df_grouped_large.loc[phylum, protein]
        x = j  # Protein position (column)
        y = i  # Phylum position (row)

        # Full circle for ≥ 50% presence
        if proportion >= 0.5:
            circle = Circle((x, y), circle_radius, edgecolor='black', facecolor=fill_color, alpha=0.8)
            ax.add_patch(circle)
        # Semi-circle for < 50% presence
        elif proportion > 0:
            half_filled_circle = Wedge((x, y), circle_radius, 0, 180, color=light_color, alpha=0.8)
            ax.add_patch(half_filled_circle)
            outline_circle = Circle((x, y), circle_radius, edgecolor='black', facecolor='none', linewidth=1)
            ax.add_patch(outline_circle)
        # Outline circle for absent proteins (0 presence)
        else:
            outline_circle = Circle((x, y), circle_radius, edgecolor='black', facecolor='none', linewidth=1)
            ax.add_patch(outline_circle)

# Step 8: Set labels with species count in parentheses
ax.set_xticks(range(len(df_grouped_large.columns)))
ax.set_xticklabels(df_grouped_large.columns, rotation=90)
ax.set_yticks(range(len(df_grouped_large.index)))
ax.set_yticklabels(phylum_labels_with_count)

# Color the y-axis labels based on domain
for i, (tick, phylum) in enumerate(zip(ax.get_yticklabels(), df_grouped_large.index)):
    domain = phylum_to_domain.get(phylum, 'Bacteria')
    tick.set_color(domain_colors.get(domain, '#CE0055'))

# Step 8: Adjust axis limits and aspect ratio
ax.set_xlim(-1, len(df_grouped_large.columns))
ax.set_ylim(-1, len(df_grouped_large.index))
ax.set_aspect('equal')

# Print the number of proteins present in each phylum
print("Number of proteins present in each phylum:")
print(df_filtered.groupby('Phylum')[protein_names_filtered].apply(lambda x: (x > 0).sum(axis=1).sum()))

# Show the plot
plt.tight_layout()
# Save the figure
plt.savefig(output_file, format = "pdf", dpi=300)
print(f"Plot saved to {output_file}")
plt.close()

