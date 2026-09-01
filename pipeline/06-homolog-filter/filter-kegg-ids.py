"""
Filter to find IDs that are NOT in KEGG list or have wrong K numbers
"""
import pandas as pd
import sys
import os


def get_expected_k_number(protein_name):
    """Get expected K number from protein_kegg.tsv"""
    try:
        protein_kegg_file = "../../input/kegg/protein_kegg.tsv"
        df = pd.read_csv(protein_kegg_file)
        match = df[df['protein'] == protein_name]
        if not match.empty:
            return match.iloc[0]['kegg']
        return None
    except Exception as e:
        print(f"Warning: Could not load protein_kegg.tsv: {e}")
        return None


def load_kegg_ids_and_numbers(kegg_file, protein_name):
    """Load KEGG IDs and determine expected K number"""
    try:
        kegg_ids = set()
        k_numbers_in_file = set()
        
        with open(kegg_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    # Only include IDs that have a second column (K number)
                    if len(parts) >= 2 and parts[1].strip():
                        kegg_ids.add(parts[0].strip())
                        k_numbers_in_file.add(parts[1].strip())
        
        # Get the expected K number from protein_kegg.tsv
        expected_k_number = get_expected_k_number(protein_name)
        
        if expected_k_number:
            print(f"Expected K number from protein_kegg.tsv: {expected_k_number}")
            return kegg_ids, {expected_k_number}
        else:
            # Fallback to all K numbers found in file
            print(f"No expected K number found, using all K numbers from file: {k_numbers_in_file}")
            return kegg_ids, k_numbers_in_file
            
    except Exception as e:
        print(f"Error loading KEGG file: {e}")
        sys.exit(1)


def filter_non_kegg_ids(table_file, kegg_ids, k_numbers, output_file):
    """Find IDs that are NOT in KEGG list OR have wrong K numbers"""
    try:
        # Read the table
        df = pd.read_csv(table_file, sep='\t')
        
        print(f"Original table size: {len(df)} rows")
        print(f"KEGG IDs: {len(kegg_ids)} IDs")
        print(f"Expected K numbers: {k_numbers}")
        
        # Filter valid IDs (not False, not missing)
        valid_df = df[(df['id'] != False) & (df['id'].notna())]
        
        # Find IDs to remove:
        # 1. IDs not in kegg_ids (not in KEGG file at all, or no K number in KEGG file)
        # 2. IDs with wrong K number (K number in table doesn't match expected K numbers)
        non_kegg_ids = []
        
        for _, row in valid_df.iterrows():
            protein_id = row['id']
            protein_kegg = row.get('kegg', False)
            
            # Skip if protein_id is False or 'False' string
            if protein_id == False or str(protein_id).upper() == 'FALSE':
                continue
            
            # Case 1: ID not in kegg_ids (either not in file or has no K number)
            if protein_id not in kegg_ids:
                non_kegg_ids.append(protein_id)
            # Case 2: ID is in kegg_ids but has wrong K number in table
            elif protein_kegg and protein_kegg != False and protein_kegg not in k_numbers:
                non_kegg_ids.append(protein_id)
                print(f"  {protein_id} has wrong K number: {protein_kegg} (expected one of {k_numbers})")
        
        # Remove duplicates and filter out False values
        non_kegg_ids = sorted(set([id for id in non_kegg_ids if id != False and str(id).upper() != 'FALSE']))
        
        print(f"Non-KEGG IDs found: {len(non_kegg_ids)} IDs")
        
        # Save IDs to text file (one per line)
        with open(output_file, 'w') as f:
            for id in non_kegg_ids:
                f.write(str(id) + '\n')
        
        print(f"Non-KEGG IDs saved to {output_file}")
        
    except Exception as e:
        print(f"Error filtering table: {e}")
        sys.exit(1)


def main(table_file, kegg_file, output_file):
    # Extract protein name from kegg_file path (e.g., TRMT112_kegg.txt -> TRMT112)
    protein_name = os.path.basename(kegg_file).replace('_kegg.txt', '')
    
    kegg_ids, k_numbers = load_kegg_ids_and_numbers(kegg_file, protein_name)
    filter_non_kegg_ids(table_file, kegg_ids, k_numbers, output_file)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter-kegg-ids.py <table_file> <kegg_file> <output_file>")
        sys.exit(1)
    
    table_file = sys.argv[1]
    kegg_file = sys.argv[2]
    output_file = sys.argv[3]
    
    main(table_file, kegg_file, output_file)
