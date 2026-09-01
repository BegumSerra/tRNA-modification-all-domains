#!/usr/bin/env python3
"""
Generate a summary table of filtering steps and sequence counts for each protein.
Shows: protein name, k01 (KEGG filter), f01/f00 (domain filter), e01/e00 (euk filter), c01/c00 (CD-HIT)
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
import pandas as pd

# Set pandas display options to show all columns
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

def get_args():
    """Parse command line arguments"""
    if len(sys.argv) < 4:
        print("Usage: python filtering_summary.py <results_dir> <blast_merged_dir> <output_tsv>")
        sys.exit(1)
    return sys.argv[1], sys.argv[2], sys.argv[3]

def count_sequences(fasta_file):
    """Count sequences in a FASTA file"""
    if not os.path.exists(fasta_file):
        return None
    try:
        return sum(1 for _ in SeqIO.parse(fasta_file, 'fasta'))
    except:
        return None

def main(results_dir, blast_merged_dir, output_file):
    
    # Get all protein names from k01 files
    proteins = []
    for f in Path(results_dir).glob("*_k01.fa"):
        protein = f.stem.replace("_k01", "")
        proteins.append(protein)
    
    proteins = sorted(set(proteins))
    
    # Collect data
    data = []
    for protein in proteins:
        row = {"Protein": protein}
        
        # Initial: from blast/merged (before KEGG filter)
        initial_file = f"{blast_merged_dir}/{protein}.fa"
        row["Initial"] = count_sequences(initial_file)
        
        # k01: KEGG filtered
        k01_file = f"{results_dir}/{protein}_k01.fa"
        row["k01"] = count_sequences(k01_file)
        
        # Determine if f01 or f00 (domain filtering)
        f01_file = f"{results_dir}/{protein}_k01_f01.fa"
        f00_file = f"{results_dir}/{protein}_k01_f00.fa"
        
        if os.path.exists(f01_file):
            row["f01"] = count_sequences(f01_file)
            row["f00"] = None
            filter_type = "f01"
        elif os.path.exists(f00_file):
            row["f00"] = count_sequences(f00_file)
            row["f01"] = None
            filter_type = "f00"
        else:
            row["f01"] = None
            row["f00"] = None
            filter_type = None
        
        # Determine if e01 or e00 (eukaryote filtering)
        if filter_type:
            e01_file = f"{results_dir}/{protein}_k01_{filter_type}_e01.fa"
            e00_file = f"{results_dir}/{protein}_k01_{filter_type}_e00.fa"
            
            if os.path.exists(e01_file):
                row["e01"] = count_sequences(e01_file)
                row["e00"] = None
                euk_type = "e01"
            elif os.path.exists(e00_file):
                row["e00"] = count_sequences(e00_file)
                row["e01"] = None
                euk_type = "e00"
            else:
                row["e01"] = None
                row["e00"] = None
                euk_type = None
            
            # c01/c00: CD-HIT clustering
            if euk_type:
                c01_file = f"{results_dir}/{protein}_k01_{filter_type}_{euk_type}_c01.fa"
                c00_file = f"{results_dir}/{protein}_k01_{filter_type}_{euk_type}_c00.fa"
                
                # Check which file exists and has sequences (>50 bytes)
                c01_count = count_sequences(c01_file) if os.path.exists(c01_file) and os.path.getsize(c01_file) > 50 else None
                c00_count = count_sequences(c00_file) if os.path.exists(c00_file) and os.path.getsize(c00_file) > 50 else None
                
                row["c01"] = c01_count
                row["c00"] = c00_count
            else:
                row["c01"] = None
                row["c00"] = None
        else:
            row["e01"] = None
            row["e00"] = None
            row["c01"] = None
            row["c00"] = None
        # Use c01 if it has sequences, otherwise use c00
        final_count = row["c01"] if row["c01"] else row["c00"]
        if row["k01"] and final_count:
            row["Total_Reduction"] = f"{100 * (1 - final_count/row['k01']):.1f}%"
        else:
            row["Total_Reduction"] = None
        
        data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Reorder columns - only include columns that exist
    columns = ["Protein", "Initial", "k01", "f01", "f00", "e01", "e00", "c01", "c00", "Total_Reduction"]
    # Filter to only existing columns
    columns_filtered = [col for col in columns if col in df.columns]
    df = df[columns_filtered]
    
    # Print summary
    print("\n" + "="*110)
    print("FILTERING SUMMARY TABLE")
    print("="*110)
    print(f"\nInitial = Original sequences from BLAST merged results (seqs/blast/merged)")
    print(f"k01 = KEGG filter (non-KEGG sequences removed)")
    print(f"f01 = Domain filter applied (sequences without domain or with extra domains removed)")
    print(f"f00 = No domain filter (all sequences kept)")
    print(f"e01 = Eukaryote filter applied (eukaryotic duplicates removed)")
    print(f"e00 = No eukaryote filter (all sequences kept)")
    print(f"c01 = CD-HIT clustering at 95% identity (sequences were reduced)")
    print(f"c00 = CD-HIT run but no clustering (all sequences kept as individual clusters)")
    print("\n" + "="*110 + "\n")
    
    # Print table
    print(df.to_string(index=False))
    
    # Save to file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\n\nTable saved to: {output_file}")
    # Print statistics
    print("\n" + "="*110)
    print("STATISTICS")
    total_proteins = len(df)
    f01_count = df['f01'].notna().sum()
    f00_count = df['f00'].notna().sum()
    e01_count = df['e01'].notna().sum()
    e00_count = df['e00'].notna().sum()
    c01_count = df['c01'].notna().sum() if 'c01' in df.columns else 0
    c00_count = df['c00'].notna().sum() if 'c00' in df.columns else 0
    
    print(f"\nTotal proteins: {total_proteins}")
    print(f"  - With domain filter (f01): {f01_count} ({100*f01_count/total_proteins:.1f}%)")
    print(f"  - Without domain filter (f00): {f00_count} ({100*f00_count/total_proteins:.1f}%)")
    print(f"  - With eukaryote filter (e01): {e01_count} ({100*e01_count/total_proteins:.1f}%)")
    print(f"  - Without eukaryote filter (e00): {e00_count} ({100*e00_count/total_proteins:.1f}%)")
    print(f"  - With CD-HIT clustering (c01): {c01_count} ({100*c01_count/total_proteins:.1f}%)")
    print(f"  - Without CD-HIT clustering (c00): {c00_count} ({100*c00_count/total_proteins:.1f}%)")
    
    print(f"\nTotal sequences:")
    if 'Initial' in df.columns and df['Initial'].notna().any():
        print(f"  - Initial (before filtering): {int(df['Initial'].sum())}")
    print(f"  - After KEGG filter (k01): {int(df['k01'].sum())}")
    # Sum c01 and c00 counts
    final_sum = (df['c01'].fillna(0).sum() + df['c00'].fillna(0).sum())
    print(f"  - After all filters (c01/c00): {int(final_sum)}")
    if 'Initial' in df.columns and df['Initial'].notna().any():
        print(f"  - Overall reduction from initial: {100 * (1 - final_sum/df['Initial'].sum()):.1f}%")
    print(f"  - Reduction from k01 to final: {100 * (1 - final_sum/df['k01'].sum()):.1f}%")
    print(f"  - Reduction from k01 to final: {100 * (1 - final_sum/df['k01'].sum()):.1f}%")
    print("\n" + "="*110 + "\n")
if __name__ == "__main__":
    if len(sys.argv) == 1:
        # Default paths for backward compatibility
        results_dir = "../../output/seqs/filtering/fasta"
        blast_merged_dir = "../../output/seqs/blast/merged"
        output_file = "../../output/table/filtering_summary.tsv"
    else:
        results_dir, blast_merged_dir, output_file = get_args()
    
    main(results_dir, blast_merged_dir, output_file)
