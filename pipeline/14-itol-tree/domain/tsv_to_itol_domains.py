#!/usr/bin/env python3
"""
Convert domain summary TSV to iTOL DATASET_DOMAINS format.

Input: TSV file with columns:
    Tree node ID, Tree node label, Protein size, Domain definitions (tab-separated)
    Domain definitions format: SHAPE|start|end|color|PF#

Output: iTOL DATASET_DOMAINS file
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


def parse_legend_file(legend_file: Path) -> dict:
    """Parse legend file to get domain info for legend block."""
    domains = {}
    with open(legend_file) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                domain, shape, color = parts[0], parts[1], parts[2]
                domains[domain] = {"shape": shape, "color": color}
    return domains


def generate_legend_block(domains: dict) -> str:
    """Generate iTOL legend block from domain info."""
    if not domains:
        return ""
    
    shapes = []
    colors = []
    labels = []
    
    for domain_name, info in sorted(domains.items()):
        shapes.append(info["shape"])
        colors.append(info["color"])
        labels.append(domain_name)
    
    legend = f"""LEGEND_TITLE,Protein Domains
LEGEND_SHAPES,{",".join(shapes)}
LEGEND_COLORS,{",".join(colors)}
LEGEND_LABELS,{",".join(labels)}
"""
    return legend


def convert_tsv_to_itol(input_tsv: Path, output_file: Path, dataset_label: str = None):
    """
    Convert domain TSV to iTOL DATASET_DOMAINS format.
    
    Args:
        input_tsv: Path to input TSV file
        output_file: Path to output iTOL file
        dataset_label: Label for the dataset (default: derived from filename)
    """
    if dataset_label is None:
        dataset_label = input_tsv.stem.replace("_domains", "")
    
    # Check for legend file
    legend_file = input_tsv.parent / f"{input_tsv.stem}_legend.tsv"
    legend_block = ""
    if legend_file.exists():
        domains = parse_legend_file(legend_file)
        legend_block = generate_legend_block(domains)
    
    # Read TSV and convert
    data_lines = []
    with open(input_tsv) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            
            node_id = parts[0]
            protein_size = parts[2]
            
            # Collect domain definitions (columns 3+)
            domain_defs = []
            for i in range(3, len(parts)):
                if parts[i].strip():
                    domain_defs.append(parts[i])
            
            if domain_defs:
                # Format: node_id,protein_size,domain1,domain2,...
                data_lines.append(f"{node_id},{protein_size},{','.join(domain_defs)}")
            else:
                # No domains, just protein backbone
                data_lines.append(f"{node_id},{protein_size}")
    
    # Write output
    with open(output_file, "w") as out:
        out.write("DATASET_DOMAINS\n")
        out.write("SEPARATOR COMMA\n")
        out.write(f"DATASET_LABEL,{dataset_label}\n")
        out.write("COLOR,#ff0000\n")
        out.write("\n")
        out.write("#optional settings\n")
        out.write("SHOW_DOMAIN_LABELS,1\n")
        out.write("BACKBONE_COLOR,#cccccc\n")
        out.write("BACKBONE_HEIGHT,6\n")
        out.write("HEIGHT_FACTOR,0.8\n")
        out.write("BORDER_WIDTH,1\n")
        out.write("BORDER_COLOR,#333333\n")
        out.write("\n")
        if legend_block:
            out.write("#legend\n")
            out.write(legend_block)
            out.write("\n")
        out.write("DATA\n")
        for line in data_lines:
            out.write(line + "\n")
    
    print(f"Converted {input_tsv} -> {output_file}")
    print(f"  {len(data_lines)} proteins")


def main():
    parser = argparse.ArgumentParser(
        description="Convert domain TSV to iTOL DATASET_DOMAINS format"
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input TSV file from create_domain_table.py"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output iTOL file (default: INPUT_itol.txt)"
    )
    parser.add_argument(
        "-l", "--label",
        help="Dataset label (default: derived from filename)"
    )
    
    args = parser.parse_args()
    
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if args.output:
        output_file = args.output
    else:
        output_file = args.input.parent / f"{args.input.stem}_itol.txt"
    
    convert_tsv_to_itol(args.input, output_file, args.label)


if __name__ == "__main__":
    main()
