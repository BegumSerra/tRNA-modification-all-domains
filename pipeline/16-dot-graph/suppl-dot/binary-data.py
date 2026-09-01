import pandas as pd
import sys

file_path = sys.argv[1]
order_file = sys.argv[2]
output_file = sys.argv[3]

df = pd.read_excel(file_path, engine="openpyxl")

# Load phylum order from file
with open(order_file, 'r') as f:
    phylum_order = [line.strip() for line in f if line.strip()]

# Define column mapping: display_name -> actual_column_name
# For proteins with multiple names (e.g., MnmE/GTPB3), use first protein's data
column_mapping = {}  # Add more mappings if needed

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
# Assume first 5 columns are: domain, phylum, order, family, species
binary_df = df.copy()

# Convert only protein data (columns after the first 5) into binary
binary_df.iloc[:, 5:] = binary_df.iloc[:, 5:].applymap(
    lambda x: 1 if pd.to_numeric(x, errors='coerce') > 0 else 0
)

# Build the final dataframe with renamed columns
result_df = pd.DataFrame()

# Keep domain column for coloring
result_df['Domain'] = binary_df['Domain']

# Define Asgard archaea orders
asgard_orders = [
    'Heimdallarchaeales',
    'Hodarchaeales',
    'Odinarchaeales',
    'Sigynarchaeales',
    'Thorarchaeales'
]

# For Archaea, use the 'order' column; for others (Eukarya, Bacteria), use 'phylum'
# For Asgard archaea orders, group them under 'Asgard'
def get_phylum_group(row):
    if row['Domain'] == 'Archaea':
        if row['Order'] in asgard_orders:
            return 'Asgard'
        return row['Order']
    return row['Phylum']
result_df['Phylum'] = binary_df.apply(get_phylum_group, axis=1)

# Add proteins in the specified order, using mapping when needed
for display_name in protein_order:
    # Get the actual column name (use mapping if exists, otherwise use display name)
    actual_col = column_mapping.get(display_name, display_name)
    
    # Add column if it exists in the data
    if actual_col in binary_df.columns:
        result_df[display_name] = binary_df[actual_col]

# Sort rows by phylum order from order.txt
# Create a mapping of phylum to order index
phylum_to_order = {phylum: i for i, phylum in enumerate(phylum_order)}
# Add order column for sorting (use large number for phyla not in order file)
result_df['_sort_order'] = result_df['Phylum'].map(lambda x: phylum_to_order.get(x, 999))
# Sort by the order
result_df = result_df.sort_values('_sort_order').drop('_sort_order', axis=1)

# Save the processed data to a CSV file with semicolon separator
result_df.to_csv(output_file, index=False, sep=';')
print(f"Binary data saved to {output_file}")
