#!/usr/bin/env python
"""
Reorder columns in the PA table Excel file
"""

import sys
import pandas as pd


def reorder_excel_columns(input_file, output_file):
    # Define the desired column order
    column_order = [
        'Domain', 'Phylum', 'Order', 'Family', 'Species',
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
    
    # Read the input file (TSV or Excel)
    print(f"Reading {input_file}...")
    if input_file.endswith('.tsv'):
        df = pd.read_csv(input_file, sep='\t')
    else:
        df = pd.read_excel(input_file, engine='openpyxl')
    
    # Get current columns
    current_columns = df.columns.tolist()
    
    # Find columns that are in the file but not in our order list
    extra_columns = [col for col in current_columns if col not in column_order]
    
    # Build final column order: specified order + any extra columns at the end
    final_column_order = []
    for col in column_order:
        if col in current_columns:
            final_column_order.append(col)
    
    # Add any extra columns at the end
    final_column_order.extend(extra_columns)
    
    # Reorder the DataFrame
    df = df[final_column_order]
    
    # Save to Excel
    print(f"Saving reordered table to {output_file}...")
    df.to_excel(output_file, index=False, engine='openpyxl')
    
    print("Done!")
    print(f"  Columns reordered: {len(final_column_order)}")
    print(f"  Rows: {len(df)}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python reorder-columns.py <input.xlsx> <output.xlsx>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    reorder_excel_columns(input_file, output_file)
