#!/usr/bin/env python3
"""
Extract protein IDs from highlighted Excel file.
Filter for sequences with kegg=False and e-value below threshold.
Can filter for archaeal only or all domains.
"""

import sys
import pandas as pd

def extract_ids(excel_file, output_file, evalue_threshold=1e-25, all_domains=False):
    """
    Extract protein IDs with kegg=False and good e-values.
    
    Args:
        excel_file: Path to highlighted Excel file
        output_file: Path to output text file with IDs
        evalue_threshold: E-value cutoff (default: 1e-25)
        all_domains: If True, include all domains; if False, only Archaea (default: False)
    """
    # Read the Excel file
    df = pd.read_excel(excel_file)
    
    # Convert e-value to float
    df['e-value_numeric'] = pd.to_numeric(df['e-value'], errors='coerce')
    
    # Handle kegg column - strip quotes and check for False
    # Excel may store as 'False, 'False', False, or "'False"
    def is_kegg_false(val):
        if pd.isna(val):
            return False
        val_str = str(val).strip().strip("'").strip('"').lower()
        return val_str == 'false'
    
    df['kegg_is_false'] = df['kegg'].apply(is_kegg_false)
    
    # Filter based on domain
    if all_domains:
        # All domains with kegg = False and e-value < threshold
        filtered = df[
            (df['kegg_is_false']) & 
            (df['e-value_numeric'] < evalue_threshold)
        ]
        domain_msg = "all domains"
    else:
        # Only archaeal species with kegg = False and e-value < threshold
        filtered = df[
            (df['domain_tax'] == 'Archaea') & 
            (df['kegg_is_false']) & 
            (df['e-value_numeric'] < evalue_threshold)
        ]
        domain_msg = "Archaea"
    
    print(f'Total {domain_msg} with kegg=False and e-value < {evalue_threshold}: {len(filtered)}')
    print('\nSample protein IDs:')
    for idx, row in filtered.head(10).iterrows():
        domain = row.get('domain_tax', 'Unknown')
        print(f"{row['id']}  {domain:10s}  {row['species']:40s}  e-value: {row['e-value']}")
        
    # Collect all IDs
    ids = filtered['id'].tolist()
    
    # Save to file
    with open(output_file, 'w') as f:
        for pid in ids:
            f.write(f"{pid}\n")
    
    print(f"\nSaved {len(ids)} IDs to {output_file}")
    return len(ids)

if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print("Usage: python find-ids.py <excel_file> <output_file> [evalue_threshold] [all_domains]")
        print("  evalue_threshold: optional, default is 1e-25")
        print("  all_domains: optional, 'all' for all domains or 'archaea' for archaea only (default: archaea)")
        sys.exit(1)
    
    excel_file = sys.argv[1]
    output_file = sys.argv[2]
    evalue_threshold = float(sys.argv[3]) if len(sys.argv) >= 4 else 1e-25
    all_domains = (sys.argv[4].lower() == 'all') if len(sys.argv) == 5 else False
    
    extract_ids(excel_file, output_file, evalue_threshold, all_domains)