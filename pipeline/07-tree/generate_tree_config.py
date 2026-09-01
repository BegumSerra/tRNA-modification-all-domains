#!/usr/bin/env python3
"""
Generate tree.yaml configuration file based on filtered fasta files.
This script scans the fasta folder and identifies the final (most filtered) version
of each protein based on the longest suffix in the filename.
"""

from pathlib import Path
from collections import defaultdict
import re

def extract_protein_name(filename):
    """Extract protein name from filename (e.g., 'ADAT2_k01_f00_e01_c01.fa' -> 'ADAT2')"""
    match = re.match(r'^([A-Za-z0-9]+)_k\d+', filename)
    if match:
        return match.group(1)
    return None

def generate_tree_config(fasta_dir, output_file):
    """
    Generate tree.yaml config file from fasta files.
    
    Args:
        fasta_dir: Path to directory containing fasta files
        output_file: Path to output YAML file
    """
    fasta_path = Path(fasta_dir)
    
    # Group files by protein name
    protein_files = defaultdict(list)
    
    for fasta_file in sorted(fasta_path.glob("*.fa")):
        protein_name = extract_protein_name(fasta_file.name)
        if protein_name:
            # Store filename without extension
            base_name = fasta_file.stem
            protein_files[protein_name].append(base_name)
    
    # For each protein, find the file with the longest name (most filtering steps)
    final_proteins = {}
    for protein, files in sorted(protein_files.items()):
        # Sort by length (longest = most filtered) and alphabetically
        final_file = max(files, key=lambda x: (len(x), x))
        final_proteins[protein] = final_file
    
    # Count statistics
    c01_count = sum(1 for f in final_proteins.values() if f.endswith('_c01'))
    c00_count = sum(1 for f in final_proteins.values() if f.endswith('_c00'))
    m01_count = sum(1 for f in final_proteins.values() if '_m01' in f)
    m02_count = sum(1 for f in final_proteins.values() if '_m02' in f)
    
    c00_proteins = [p for p, f in final_proteins.items() if f.endswith('_c00')]
    m02_proteins = [p for p, f in final_proteins.items() if '_m02' in f]
    
    # Generate YAML content
    yaml_lines = [
        "# Configuration for phylogenetic tree construction",
        "# Lists all proteins with their final filtered sequences (c00 or c01)",
        "PROTEIN:"
    ]
    
    # Add protein list
    protein_list = sorted(final_proteins.values())
    yaml_lines.append("  [" + '"{}"'.format(protein_list[0]) + ",")
    for i in range(1, len(protein_list)):
        yaml_lines.append(f'  "{protein_list[i]}"' + ("," if i < len(protein_list) - 1 else "]"))
    
    # Add mapping section
    yaml_lines.extend([
        "",
        "# Mapping from protein name to filter string (for file lookup)",
        "FILTER_TO_PROTEIN:"
    ])
    
    for protein in sorted(final_proteins.keys()):
        yaml_lines.append(f"  {protein}: {final_proteins[protein]}")
    
    # Add summary statistics
    yaml_lines.extend([
        "",
        "# Summary statistics",
        f"# Total proteins: {len(final_proteins)}",
        f"# With CD-HIT clustering (c01): {c01_count} proteins",
        f"# Without CD-HIT clustering (c00): {c00_count} proteins" + (f" ({', '.join(c00_proteins)})" if c00_proteins else ""),
        f"# With manual filtering (m01): {m01_count} proteins",
        f"# With manual filtering (m02): {m02_count} proteins" + (f" ({', '.join(m02_proteins)})" if m02_proteins else ""),
        "",
        "# Parameters",
        "uf: 1000",
        ""
    ])
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(yaml_lines))
    
    print(f"✓ Generated {output_file}")
    print(f"  Total proteins: {len(final_proteins)}")
    print(f"  With c01 (clustered): {c01_count}")
    print(f"  With c00 (not clustered): {c00_count}")
    if m02_count > 0:
        print(f"  With manual filtering (m02): {m02_count}")

if __name__ == "__main__":
    # Default paths (relative to script location)
    script_dir = Path(__file__).parent
    fasta_dir = script_dir / "../../output/seqs/filtering/fasta"
    output_file = script_dir / "tree.yaml"
    
    # Allow command-line override
    import sys
    if len(sys.argv) > 1:
        fasta_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_file = Path(sys.argv[2])
    
    generate_tree_config(fasta_dir, output_file)
