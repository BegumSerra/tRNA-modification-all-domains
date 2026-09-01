#!/usr/bin/env python3
"""
Color only presence (value=1) cells in the grouped presence/absence Excel table.
"""

import sys
from openpyxl import load_workbook
from openpyxl.styles import PatternFill


def color_presence_only(input_file, output_file):
    """
    Load the grouped Excel file and color only cells with value 1.
    Cells with 0/FALSE get white background.
    
    Args:
        input_file: Path to input Excel file (protein_presence_grouped.xlsx)
        output_file: Path to output Excel file
    """
    # Load the workbook
    wb = load_workbook(input_file)
    ws = wb.active
    
    print(f"Processing {input_file}")
    print(f"Rows: {ws.max_row}, Columns: {ws.max_column}")
    
    # Get header row to identify protein columns
    headers = [ws.cell(row=1, column=col).value for col in range(1, ws.max_column + 1)]
    taxonomy_cols = ['Domain', 'Phylum', 'Order', 'Family', 'Species']
    protein_col_indices = [i + 1 for i, h in enumerate(headers) if h not in taxonomy_cols]
    
    print(f"Protein columns: {len(protein_col_indices)}")
    
    white_fill = PatternFill(start_color="FFFFFF", end_color="FFFFFF", fill_type="solid")
    
    cells_modified = 0
    
    # Process each row (skip header)
    for row_idx in range(2, ws.max_row + 1):
        for col_idx in protein_col_indices:
            cell = ws.cell(row=row_idx, column=col_idx)
            value = cell.value
            
            # Get current fill color
            current_fill = cell.fill
            current_color = None
            if current_fill and current_fill.fgColor and current_fill.fgColor.rgb:
                current_color = current_fill.fgColor.rgb
            
            # Skip if it's a separator row (black) or no color
            if current_color in ['00000000', '000000', None]:
                continue
            
            # Check if value is presence (1 or TRUE)
            is_presence = False
            try:
                if value is not None:
                    if isinstance(value, (int, float)):
                        is_presence = value == 1 or value >= 1
                    elif isinstance(value, str):
                        is_presence = value.strip().upper() == 'TRUE' or value.strip() == '1'
                        # Also check for numeric strings
                        try:
                            is_presence = float(value) >= 1
                        except ValueError:
                            pass
            except (ValueError, TypeError):
                pass
            
            # If not presence (i.e., 0, FALSE, or empty), set to white
            if not is_presence:
                cell.fill = white_fill
                cells_modified += 1
    
    # Save the modified workbook
    wb.save(output_file)
    print(f"Cells set to white (non-presence): {cells_modified}")
    print(f"Output saved to: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python color-presence-only.py <input_xlsx> <output_xlsx>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    color_presence_only(input_file, output_file)
