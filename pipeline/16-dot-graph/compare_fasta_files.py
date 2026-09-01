#!/usr/bin/env python3
"""
Compare FASTA files between new and old folders to find missing or extra sequences.

New files: 
  - Homologs: /home/begum/PycharmProjects/PythonProject/Jan_2026_tRNA/output/seqs/homologs
  - Combined: /home/begum/PycharmProjects/PythonProject/Jan_2026_tRNA/output/seqs/combined_tree
Old files: /home/begum/Projects/EBA/output/seqs/filtering

New file naming: 
- For most proteins: use r01 or r00 suffix (e.g., ADAT2_k01_f01_e01_c01_nogaps_r01.fa)
- For Tuc1 and TRM9B: use a01 suffix (e.g., Tuc1_k01_f01_e00_c01_nogaps_r01_a01.fa)
- For combined proteins: use m01_r01 suffix (e.g., TadA-ADAT2_m01_r01.fa)

Old file naming priority: m02, m01, e02, c01, e01 (use most filtered version available)
"""

import os
import glob
import re
from collections import defaultdict

# Directories
NEW_DIR = "/home/begum/PycharmProjects/PythonProject/Jan_2026_tRNA/output/seqs/homologs"
NEW_COMBINED_DIR = "/home/begum/PycharmProjects/PythonProject/Jan_2026_tRNA/output/seqs/combined_tree"
OLD_DIR = "/home/begum/Projects/EBA/output/seqs/filtering"

# Special proteins that use a01 suffix instead of r01/r00
SPECIAL_PROTEINS = ['Tuc1', 'TRM9B']

# Combined proteins (name mapping: new_name -> old_name patterns)
COMBINED_PROTEINS = {
    'TadA-ADAT2': ['Tada-ADAT2', 'TadA-ADAT2'],
    'MnmA-MTU1': ['MnmA-MTU1'],
    'Tgt-QTRT1': ['Tgt-QTRT1', 'TgtA-Tgt-QTRT1']
}

def parse_fasta(filepath):
    """Parse a FASTA file and return a dict of {sequence_id: sequence}."""
    sequences = {}
    current_id = None
    current_seq = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    # Get just the ID part (first word after >)
                    current_id = line[1:].split()[0]
                    current_seq = []
                elif line:
                    current_seq.append(line)
            
            if current_id:
                sequences[current_id] = ''.join(current_seq)
    except FileNotFoundError:
        return None
    
    return sequences

def get_protein_name(filename):
    """Extract protein name from filename."""
    # Handle combined protein names like TadA-ADAT2
    match = re.match(r'^([A-Za-z0-9]+-[A-Za-z0-9]+)', filename)
    if match:
        return match.group(1)
    # Match single protein name at the start (before _k or _f or _nogaps)
    match = re.match(r'^([A-Za-z0-9]+)', filename)
    if match:
        return match.group(1)
    return None


def find_new_combined_file(protein_name):
    """Find the appropriate new combined FASTA file for a protein."""
    if not os.path.exists(NEW_COMBINED_DIR):
        return None
    
    files = os.listdir(NEW_COMBINED_DIR)
    
    # Look for m01_r01 suffix first, then m01, then base
    for suffix in ['_m01_r01.fa', '_m01.fa', '.fa']:
        target = f"{protein_name}{suffix}"
        if target in files:
            return os.path.join(NEW_COMBINED_DIR, target)
    
    return None


def find_old_combined_file(protein_name):
    """Find the appropriate old combined FASTA file for a protein."""
    files = os.listdir(OLD_DIR)
    
    # Get possible old names for this combined protein
    old_names = COMBINED_PROTEINS.get(protein_name, [protein_name])
    
    # Score each file based on which suffixes it contains (higher = more filtered)
    def score_file(filename):
        score = 0
        fname = filename.lower()
        if '_m03' in fname:
            score += 200
        if '_m02' in fname:
            score += 100
        if '_m01' in fname:
            score += 50
        if '_e02' in fname:
            score += 25
        if '_e01' in fname:
            score += 10
        return score
    
    # Find matching files for any of the old names
    matching_files = []
    for old_name in old_names:
        for f in files:
            if f.endswith('.fa'):
                fname_lower = f.lower()
                oname_lower = old_name.lower()
                if fname_lower.startswith(oname_lower + '_') or fname_lower == oname_lower + '.fa':
                    matching_files.append(f)
    
    if not matching_files:
        return None
    
    # Sort by score (highest first)
    matching_files.sort(key=lambda f: (score_file(f), len(f)), reverse=True)
    
    return os.path.join(OLD_DIR, matching_files[0])


def find_new_file(protein_name):
    """Find the appropriate new FASTA file for a protein."""
    files = os.listdir(NEW_DIR)
    
    # For special proteins, look for a01 suffix
    if protein_name in SPECIAL_PROTEINS:
        pattern = f"{protein_name}_*_a01.fa"
        matches = [f for f in files if re.match(f"^{protein_name}_.*_a01\\.fa$", f)]
        if matches:
            return os.path.join(NEW_DIR, matches[0])
    
    # For other proteins, look for r01 first, then r00
    for suffix in ['_r01.fa', '_r00.fa']:
        matches = [f for f in files if f.startswith(protein_name + '_') and f.endswith(suffix)]
        if matches:
            # Prefer the longest match (most specific)
            matches.sort(key=len, reverse=True)
            return os.path.join(NEW_DIR, matches[0])
    
    return None

def find_old_file(protein_name):
    """Find the appropriate old FASTA file for a protein using priority: m02, m01, e02, c01, e01."""
    files = os.listdir(OLD_DIR)
    
    # Priority order for suffixes (most filtered first)
    priority_suffixes = ['_m02.fa', '_m01.fa', '_e02.fa', '_c01.fa', '_e01.fa']
    
    # Also handle combined suffixes like _e01_c01_m02.fa
    protein_files = [f for f in files if f.startswith(protein_name + '_') or f.startswith(protein_name + '.') or f == f"{protein_name}_nogaps.fa"]
    protein_files = [f for f in files if f.lower().startswith(protein_name.lower() + '_') or 
                     f.lower().startswith(protein_name.lower() + '.') or
                     f.lower() == f"{protein_name.lower()}_nogaps.fa"]
    
    # Try to find files with case-insensitive match
    protein_files = []
    for f in files:
        fname_lower = f.lower()
        pname_lower = protein_name.lower()
        if fname_lower.startswith(pname_lower + '_') and f.endswith('.fa'):
            protein_files.append(f)
    
    if not protein_files:
        return None
    
    # Score each file based on which suffixes it contains (higher = more filtered)
    def score_file(filename):
        score = 0
        fname = filename.lower()
        if '_m02' in fname:
            score += 100
        if '_m01' in fname:
            score += 50
        if '_e02' in fname:
            score += 25
        if '_c01' in fname:
            score += 10
        if '_e01' in fname:
            score += 5
        return score
    
    # Sort by score (highest first) then by length (longest first for ties)
    protein_files.sort(key=lambda f: (score_file(f), len(f)), reverse=True)
    
    return os.path.join(OLD_DIR, protein_files[0])

def get_all_proteins():
    """Get list of all protein names from both directories."""
    new_files = os.listdir(NEW_DIR)
    old_files = os.listdir(OLD_DIR)
    
    proteins = set()
    
    # Extract protein names from new files
    for f in new_files:
        if f.endswith('.fa'):
            name = get_protein_name(f)
            if name:
                proteins.add(name)
    
    # Extract protein names from old files (case-insensitive matching)
    for f in old_files:
        if f.endswith('.fa'):
            name = get_protein_name(f)
            if name:
                proteins.add(name)
    
    return sorted(proteins)

def compare_sequences(new_seqs, old_seqs):
    """Compare two sets of sequences and return differences."""
    new_ids = set(new_seqs.keys())
    old_ids = set(old_seqs.keys())
    
    only_in_new = new_ids - old_ids
    only_in_old = old_ids - new_ids
    in_both = new_ids & old_ids
    
    # Check for sequence differences in shared IDs
    seq_differences = []
    for seq_id in in_both:
        if new_seqs[seq_id] != old_seqs[seq_id]:
            seq_differences.append(seq_id)
    
    return {
        'only_in_new': only_in_new,
        'only_in_old': only_in_old,
        'in_both': in_both,
        'seq_differences': seq_differences
    }

def main():
    proteins = get_all_proteins()
    
    print("=" * 80)
    print("FASTA FILE COMPARISON: NEW vs OLD (SINGLE PROTEINS)")
    print("=" * 80)
    print(f"\nNew directory: {NEW_DIR}")
    print(f"Old directory: {OLD_DIR}")
    print(f"\nTotal proteins found: {len(proteins)}")
    print(f"Special proteins (using a01): {SPECIAL_PROTEINS}")
    print("\n" + "=" * 80)
    
    results = {}
    missing_new = []
    missing_old = []
    
    for protein in proteins:
        new_file = find_new_file(protein)
        old_file = find_old_file(protein)
        
        if not new_file:
            missing_new.append(protein)
            continue
        
        if not old_file:
            missing_old.append(protein)
            continue
        
        new_seqs = parse_fasta(new_file)
        old_seqs = parse_fasta(old_file)
        
        if new_seqs is None:
            missing_new.append(protein)
            continue
        
        if old_seqs is None:
            missing_old.append(protein)
            continue
        
        comparison = compare_sequences(new_seqs, old_seqs)
        results[protein] = {
            'new_file': os.path.basename(new_file),
            'old_file': os.path.basename(old_file),
            'new_count': len(new_seqs),
            'old_count': len(old_seqs),
            'comparison': comparison
        }
    
    # Print summary for single proteins
    print("\n" + "=" * 80)
    print("MISSING FILES (SINGLE PROTEINS)")
    print("=" * 80)
    
    if missing_new:
        print(f"\nProteins without NEW file ({len(missing_new)}):")
        for p in missing_new:
            print(f"  - {p}")
    else:
        print("\nAll proteins have NEW files ✓")
    
    if missing_old:
        print(f"\nProteins without OLD file ({len(missing_old)}):")
        for p in missing_old:
            print(f"  - {p}")
    else:
        print("\nAll proteins have OLD files ✓")
    
    # Print detailed comparison for single proteins
    print("\n" + "=" * 80)
    print("SEQUENCE COMPARISON (SINGLE PROTEINS)")
    print("=" * 80)
    
    proteins_with_differences = []
    proteins_identical = []
    
    for protein, data in sorted(results.items()):
        comp = data['comparison']
        has_diff = (len(comp['only_in_new']) > 0 or 
                    len(comp['only_in_old']) > 0 or 
                    len(comp['seq_differences']) > 0)
        
        if has_diff:
            proteins_with_differences.append((protein, data))
        else:
            proteins_identical.append(protein)
    
    print(f"\n📊 Summary:")
    print(f"  - Proteins with IDENTICAL sequences: {len(proteins_identical)}")
    print(f"  - Proteins with DIFFERENCES: {len(proteins_with_differences)}")
    
    if proteins_identical:
        print(f"\n✓ Identical proteins ({len(proteins_identical)}):")
        for p in proteins_identical:
            data = results[p]
            print(f"  {p}: {data['new_count']} sequences (new: {data['new_file']}, old: {data['old_file']})")
    
    if proteins_with_differences:
        print(f"\n⚠ Proteins with differences ({len(proteins_with_differences)}):")
        print("-" * 80)
        
        for protein, data in proteins_with_differences:
            comp = data['comparison']
            print(f"\n{protein}:")
            print(f"  New file: {data['new_file']} ({data['new_count']} seqs)")
            print(f"  Old file: {data['old_file']} ({data['old_count']} seqs)")
            
            if comp['only_in_new']:
                print(f"  🆕 Only in NEW ({len(comp['only_in_new'])}):")
                for seq_id in sorted(comp['only_in_new'])[:10]:
                    print(f"      + {seq_id}")
                if len(comp['only_in_new']) > 10:
                    print(f"      ... and {len(comp['only_in_new']) - 10} more")
            
            if comp['only_in_old']:
                print(f"  🔴 Only in OLD ({len(comp['only_in_old'])}):")
                for seq_id in sorted(comp['only_in_old'])[:10]:
                    print(f"      - {seq_id}")
                if len(comp['only_in_old']) > 10:
                    print(f"      ... and {len(comp['only_in_old']) - 10} more")
            
            if comp['seq_differences']:
                print(f"  🔄 Different sequences ({len(comp['seq_differences'])}):")
                for seq_id in sorted(comp['seq_differences'])[:5]:
                    print(f"      ~ {seq_id}")
                if len(comp['seq_differences']) > 5:
                    print(f"      ... and {len(comp['seq_differences']) - 5} more")
    
    # ========================================================================
    # COMBINED PROTEINS COMPARISON
    # ========================================================================
    print("\n" + "=" * 80)
    print("FASTA FILE COMPARISON: NEW vs OLD (COMBINED PROTEINS)")
    print("=" * 80)
    print(f"\nNew combined directory: {NEW_COMBINED_DIR}")
    print(f"Old directory: {OLD_DIR}")
    print(f"\nCombined proteins to compare: {list(COMBINED_PROTEINS.keys())}")
    print("\n" + "=" * 80)
    
    combined_results = {}
    combined_missing_new = []
    combined_missing_old = []
    
    for protein_name in COMBINED_PROTEINS.keys():
        new_file = find_new_combined_file(protein_name)
        old_file = find_old_combined_file(protein_name)
        
        if not new_file:
            combined_missing_new.append(protein_name)
            continue
        
        if not old_file:
            combined_missing_old.append(protein_name)
            continue
        
        new_seqs = parse_fasta(new_file)
        old_seqs = parse_fasta(old_file)
        
        if new_seqs is None:
            combined_missing_new.append(protein_name)
            continue
        
        if old_seqs is None:
            combined_missing_old.append(protein_name)
            continue
        
        comparison = compare_sequences(new_seqs, old_seqs)
        combined_results[protein_name] = {
            'new_file': os.path.basename(new_file),
            'old_file': os.path.basename(old_file),
            'new_count': len(new_seqs),
            'old_count': len(old_seqs),
            'comparison': comparison
        }
    
    # Print summary for combined proteins
    print("\n" + "=" * 80)
    print("MISSING FILES (COMBINED PROTEINS)")
    print("=" * 80)
    
    if combined_missing_new:
        print(f"\nCombined proteins without NEW file ({len(combined_missing_new)}):")
        for p in combined_missing_new:
            print(f"  - {p}")
    else:
        print("\nAll combined proteins have NEW files ✓")
    
    if combined_missing_old:
        print(f"\nCombined proteins without OLD file ({len(combined_missing_old)}):")
        for p in combined_missing_old:
            print(f"  - {p}")
    else:
        print("\nAll combined proteins have OLD files ✓")
    
    # Print detailed comparison for combined proteins
    print("\n" + "=" * 80)
    print("SEQUENCE COMPARISON (COMBINED PROTEINS)")
    print("=" * 80)
    
    combined_with_differences = []
    combined_identical = []
    
    for protein, data in sorted(combined_results.items()):
        comp = data['comparison']
        has_diff = (len(comp['only_in_new']) > 0 or 
                    len(comp['only_in_old']) > 0 or 
                    len(comp['seq_differences']) > 0)
        
        if has_diff:
            combined_with_differences.append((protein, data))
        else:
            combined_identical.append(protein)
    
    print(f"\n📊 Summary (Combined Proteins):")
    print(f"  - Combined proteins with IDENTICAL sequences: {len(combined_identical)}")
    print(f"  - Combined proteins with DIFFERENCES: {len(combined_with_differences)}")
    
    if combined_identical:
        print(f"\n✓ Identical combined proteins ({len(combined_identical)}):")
        for p in combined_identical:
            data = combined_results[p]
            print(f"  {p}: {data['new_count']} sequences (new: {data['new_file']}, old: {data['old_file']})")
    
    if combined_with_differences:
        print(f"\n⚠ Combined proteins with differences ({len(combined_with_differences)}):")
        print("-" * 80)
        
        for protein, data in combined_with_differences:
            comp = data['comparison']
            print(f"\n{protein}:")
            print(f"  New file: {data['new_file']} ({data['new_count']} seqs)")
            print(f"  Old file: {data['old_file']} ({data['old_count']} seqs)")
            
            if comp['only_in_new']:
                print(f"  🆕 Only in NEW ({len(comp['only_in_new'])}):")
                for seq_id in sorted(comp['only_in_new'])[:10]:
                    print(f"      + {seq_id}")
                if len(comp['only_in_new']) > 10:
                    print(f"      ... and {len(comp['only_in_new']) - 10} more")
            
            if comp['only_in_old']:
                print(f"  🔴 Only in OLD ({len(comp['only_in_old'])}):")
                for seq_id in sorted(comp['only_in_old'])[:10]:
                    print(f"      - {seq_id}")
                if len(comp['only_in_old']) > 10:
                    print(f"      ... and {len(comp['only_in_old']) - 10} more")
            
            if comp['seq_differences']:
                print(f"  🔄 Different sequences ({len(comp['seq_differences'])}):")
                for seq_id in sorted(comp['seq_differences'])[:5]:
                    print(f"      ~ {seq_id}")
                if len(comp['seq_differences']) > 5:
                    print(f"      ... and {len(comp['seq_differences']) - 5} more")
    
    print("\n" + "=" * 80)
    print("COMPARISON COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()