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
df_large = pd.read_csv(input_file, header=0, sep=";")

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

# Define eukaryotic supergroups for grouping
tsar_phyla = ['Alveolata', 'Stramenopiles', 'Rhizaria', 'Telonemia']
haptista_phyla = ['Haptophyta', 'Centroplasthelida']
cryptista_phyla = ['Cryptophyceae', 'Palpitomonas']
archeaplastida_phyla = ['Chloroplastida', 'Rhodophyta', 'Glaucophyta', 'Rhodelphis', 'Ancoracysta']
amorphea_phyla = ['Apusomonadida', 'Breviatea', 'Opisthokonta', 'Amoebozoa']
crums_phyla = ['Collodictyonidae', 'Rigifilida', 'Mantamonas']
discoba_phyla = ['Euglenozoa', 'Heterolobosea', 'Jakobida']
metamonada_phyla = ['Fornicata', 'Parabasalia', 'Anaeramoebidae', 'Barthelona', 'Preaxostyla']
malawimonada_phyla = ['Malawimonadida']
ancyromonadida_phyla = ['Ancyromonadida']
hemimastigophora_phyla = ['Hemimastigophora']

def get_eukaryote_group(phylum):
    """Map a phylum to its eukaryotic supergroup."""
    if phylum in tsar_phyla:
        return 'TSAR'
    elif phylum in haptista_phyla:
        return 'Haptista'
    elif phylum in cryptista_phyla:
        return 'Cryptista'
    elif phylum in archeaplastida_phyla:
        return 'Archaeplastida'
    elif phylum in amorphea_phyla:
        return 'Amorphea'
    elif phylum in crums_phyla:
        return 'CRuMs'
    elif phylum in discoba_phyla:
        return 'Discoba'
    elif phylum in metamonada_phyla:
        return 'Metamonada'
    elif phylum in malawimonada_phyla:
        return 'Malawimonada'
    elif phylum in ancyromonadida_phyla:
        return 'Ancyromonadida'
    elif phylum in hemimastigophora_phyla:
        return 'Hemimastigophora'
    else:
        return phylum  # Return original if not in any group

# Create a mapping from phylum to domain
phylum_to_domain = df_large.groupby('Phylum')['Domain'].first().to_dict()

# Ensure the required columns exist
if 'Phylum' not in df_large.columns:
    raise ValueError("The CSV file must have a 'Phylum' column.")
if 'Domain' not in df_large.columns:
    raise ValueError("The CSV file must have a 'Domain' column.")

# Handle missing values in the 'Phylum' column: replace NaN with 'Unknown'
df_large['Phylum'] = df_large['Phylum'].fillna('Unknown')

# Apply eukaryotic grouping - replace phyla with their supergroup for Eukarya
df_large['Phylum_Grouped'] = df_large.apply(
    lambda row: get_eukaryote_group(row['Phylum']) if row['Domain'] == 'Eukarya' else row['Phylum'],
    axis=1
)

# Extract protein names (all columns except 'Domain', 'Phylum', and 'Phylum_Grouped')
protein_names = [col for col in df_large.columns if col not in ['Domain', 'Phylum', 'Phylum_Grouped']]

# Step 1: Create a list of proteins with a non-zero presence across all species
valid_proteins = [protein for protein in protein_names if df_large[protein].sum() > 0]

# Step 2: Filter the dataframe to only include the valid proteins
df_filtered = df_large[['Domain', 'Phylum', 'Phylum_Grouped'] + valid_proteins]

# Update the protein names list after filtering
protein_names_filtered = [col for col in df_filtered.columns if col not in ['Domain', 'Phylum', 'Phylum_Grouped']]

# Step 3: Count species per grouped phylum
species_per_phylum = df_filtered.groupby("Phylum_Grouped").size()

# Create mapping from grouped phylum to domain (use first occurrence)
grouped_phylum_to_domain = df_filtered.groupby('Phylum_Grouped')['Domain'].first().to_dict()

# Step 4: Group by Phylum_Grouped and calculate the proportion of presence for each protein
df_grouped_large = df_filtered.groupby("Phylum_Grouped")[protein_names_filtered].apply(lambda x: (x > 0).sum() / len(x))

# Step 5: Preserve the phylum order from order.txt
# Filter to only include phyla that exist in the data
unique_phyla_order = [p for p in phylum_order_from_file if p in df_grouped_large.index]
df_grouped_large = df_grouped_large.loc[unique_phyla_order]

# Step 6: Create phylum labels with species count in parentheses
phylum_labels_with_count = [f"{phylum} ({species_per_phylum[phylum]})" for phylum in df_grouped_large.index]

# Step 7: Create the figure and axis
fig, ax = plt.subplots(figsize=(15, 10))

# Define circle and arc sizes
circle_radius = 0.3

# Step 8: Iterate over phyla (rows) and proteins (columns)
for i, phylum in enumerate(df_grouped_large.index):
    # Get domain for this grouped phylum
    domain = grouped_phylum_to_domain.get(phylum, 'Bacteria')
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

# Step 9: Set labels with species count in parentheses
ax.set_xticks(range(len(df_grouped_large.columns)))
ax.set_xticklabels(df_grouped_large.columns, rotation=90)
ax.set_yticks(range(len(df_grouped_large.index)))
ax.set_yticklabels(phylum_labels_with_count)

# Color the y-axis labels based on domain
for i, (tick, phylum) in enumerate(zip(ax.get_yticklabels(), df_grouped_large.index)):
    domain = grouped_phylum_to_domain.get(phylum, 'Bacteria')
    tick.set_color(domain_colors.get(domain, '#CE0055'))

# Step 10: Adjust axis limits and aspect ratio
ax.set_xlim(-1, len(df_grouped_large.columns))
ax.set_ylim(-1, len(df_grouped_large.index))
ax.set_aspect('equal')

# Print the number of proteins present in each phylum
print("Number of proteins present in each phylum:")
print(df_filtered.groupby('Phylum_Grouped')[protein_names_filtered].apply(lambda x: (x > 0).sum(axis=1).sum()))

# Show the plot
plt.tight_layout()
# Save the figure
plt.savefig(output_file, format="pdf", dpi=300)
print(f"Plot saved to {output_file}")
plt.close()