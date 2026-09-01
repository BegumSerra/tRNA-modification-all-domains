import pandas as pd
import sys

file_path = sys.argv[1]
order_file = sys.argv[2]
output_file = sys.argv[3]

df = pd.read_excel(file_path, engine="openpyxl")

# Load phylum order from file
with open(order_file, 'r') as f:
    phylum_order = [line.strip() for line in f if line.strip()]

print(f"Original data shape: {df.shape}")
print(f"Original columns: {list(df.columns)}")
print(f"First few phylum values: {df['Phylum'].head()}")

# Define column mapping: display_name -> actual_column_name
column_mapping = {}  # Remove the problematic mappings since columns are already correctly named

# Define the expected column order (for display)
protein_order = [
    'MnmE/GTPB3', 'MnmG/MTO1', 'MnmC', 'MnmA-MTU1', 'IscS/NFS1',
    'TusA', 'TusC', 'TusD', 'TusE', 'SelD',
    'MnmH', 'CmoA', 'CmoB', 'CmoM',
    'ELP1', 'ELP2', 'ELP3', 'ELP4', 'ELP5', 'ELP6',
    'ALKB8', 'TRM9B', 'Trm9', 'TRMT112',
    'Tuc1/CTU1', 'CTU2', 'TrhP', 'TrhO',
    'PUS1', 'FolE', 'QueD', 'QueE', 'QueC', 'QueF', 'QueA', 'QueG',
    'Tgt-QTRT1', 'QTRT2',
    'ALKB1', 'NSUN3', 'NSUN2', 'TRM7', 'WDR6',
    'TiaS', 'TilS', 'TmcA', 'TrmL',
    'TadA-ADAT2', 'ADAT3', 'YadB'
]

# Process the data
binary_df = df.copy()

# Identify taxonomy columns (first 5 columns)
taxonomy_columns = ['Domain', 'Phylum', 'Order', 'Family', 'Species']

# Get protein column names (everything after taxonomy columns)
protein_columns = [col for col in binary_df.columns if col not in taxonomy_columns]

print(f"Taxonomy columns: {taxonomy_columns}")
print(f"Found {len(protein_columns)} protein columns")
print(f"Sample protein columns: {protein_columns[:10]}")

# Convert ONLY protein columns to binary (preserve taxonomy columns)
for col in protein_columns:
    if col in binary_df.columns:
        binary_df[col] = binary_df[col].apply(
            lambda x: 1 if (pd.notna(x) and str(x).upper() not in ['FALSE', 'N/A', ''] and pd.to_numeric(x, errors='coerce') != 0) else 0
        )

print(f"After binary conversion - Phylum values: {binary_df['Phylum'].head()}")
print(f"Unique phyla: {binary_df['Phylum'].unique()[:10]}")

# Count proteins by phylum
print("Number of species in each phylum:")
phylum_counts = binary_df['Phylum'].value_counts()
print(phylum_counts)

# Order phylogenetically if phylum order is provided
if phylum_order:
    print(f"Applying phylum order from file: {len(phylum_order)} phylums specified")
    # Create a mapping of phylum to order index
    phylum_order_dict = {phylum: idx for idx, phylum in enumerate(phylum_order)}
    
    # Add order column for sorting
    binary_df['phylum_order'] = binary_df['Phylum'].map(phylum_order_dict)
    
    # Sort by phylum order (unknown phylums will have NaN and go to the end)
    binary_df = binary_df.sort_values(['phylum_order', 'Phylum', 'Species'], na_last=True)
    
    # Remove the temporary order column
    binary_df = binary_df.drop('phylum_order', axis=1)
    
    print("Data ordered by phylum")
else:
    # Sort alphabetically by phylum and species
    binary_df = binary_df.sort_values(['Phylum', 'Species'])
    print("Data sorted alphabetically")

# Create result dataframe with taxonomy columns first
result_df = pd.DataFrame()

# Add taxonomy columns
for col in taxonomy_columns:
    if col in binary_df.columns:
        result_df[col] = binary_df[col]

# Add proteins in the specified order, using mapping when needed
proteins_added = 0
proteins_missing = []

for display_name in protein_order:
    # Get the actual column name (use mapping if exists, otherwise use display name)
    actual_col = column_mapping.get(display_name, display_name)
    
    # Add column if it exists in the data
    if actual_col in binary_df.columns:
        result_df[display_name] = binary_df[actual_col]
        proteins_added += 1
    else:
        proteins_missing.append(display_name)
        print(f"Warning: Column '{actual_col}' not found in data (display name: '{display_name}')")

# Add any remaining protein columns that weren't in the specified order
remaining_proteins = []
for col in protein_columns:
    if col not in [column_mapping.get(p, p) for p in protein_order]:
        if col in binary_df.columns:
            result_df[col] = binary_df[col]
            remaining_proteins.append(col)
            proteins_added += 1

print(f"Added {proteins_added} protein columns to result")
if proteins_missing:
    print(f"Missing proteins from order: {proteins_missing}")
if remaining_proteins:
    print(f"Additional proteins not in order: {remaining_proteins}")

# Final data validation
print(f"Final data shape: {result_df.shape}")
print(f"Final columns: {list(result_df.columns)}")

# Check for any remaining non-binary values in protein columns
final_protein_cols = [col for col in result_df.columns if col not in taxonomy_columns]
for col in final_protein_cols[:5]:  # Check first 5 protein columns
    unique_vals = result_df[col].unique()
    if len(unique_vals) > 2 or not all(val in [0, 1] for val in unique_vals):
        print(f"Warning: Column '{col}' has non-binary values: {unique_vals}")

# Save the processed data
result_df.to_csv(output_file, sep='\t', index=False)
print(f"Binary data saved to: {output_file}")

# Print summary statistics
total_proteins = len(final_protein_cols)
total_species = len(result_df)
print(f"\nSummary:")
print(f"Total species: {total_species}")
print(f"Total proteins: {total_proteins}")

# Calculate presence statistics by phylum
print(f"\nPresence statistics by phylum:")
for phylum in result_df['Phylum'].unique()[:10]:  # Show first 10 phylums
    phylum_data = result_df[result_df['Phylum'] == phylum]
    protein_presence = phylum_data[final_protein_cols].sum(axis=1).mean()
    print(f"  {phylum}: {len(phylum_data)} species, avg {protein_presence:.1f} proteins per species")