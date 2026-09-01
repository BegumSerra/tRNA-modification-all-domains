#!/usr/bin/env python3
"""
Create domain summary table from InterProScan TSV output.
Extracts Pfam domains and formats them for iTOL visualization.

Reads a FASTA file to filter proteins and uses the FASTA IDs for output.
Automatically assigns different colors and shapes to each unique domain.

Output format:
Tree node ID    Tree node label    Protein size    Domain definitions
(Domain definitions: SHAPE|start|end|color|PF#)
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


# Available shapes for iTOL domain visualization
SHAPES = ["RE", "HH", "HV", "EL", "DI", "TR", "TL", "PL", "PR", "PU", "PD", "OC", "GP"]
# RE = rectangle, HH = horizontal hexagon, HV = vertical hexagon
# EL = ellipse, DI = diamond, TR = triangle right, TL = triangle left
# PL = pentagon left, PR = pentagon right, PU = pentagon up, PD = pentagon down
# OC = octagon, GP = rounded rectangle

# Color palette - pale/pastel colors for domains
COLORS = [
    "#FFB3BA",  # pale pink
    "#BAFFC9",  # pale green
    "#BAE1FF",  # pale blue
    "#FFDFBA",  # pale orange
    "#E0BBE4",  # pale purple
    "#C9FFE5",  # pale mint
    "#FFC6FF",  # pale magenta
    "#FFFFBA",  # pale yellow
    "#FFD1DC",  # pale rose
    "#B5EAD7",  # pale teal
    "#C7CEEA",  # pale lavender
    "#E2CFC4",  # pale taupe
    "#F0E6EF",  # pale mauve
    "#D4A5A5",  # pale coral
    "#A8E6CF",  # pale seafoam
    "#DCEDC1",  # pale lime
    "#FFEEAD",  # pale cream
    "#D5AAFF",  # pale violet
    "#B8D4E3",  # pale sky
    "#F9D5E5",  # pale blush
]


def parse_fasta_ids(fasta_file: Path) -> set:
    """
    Parse FASTA file and extract sequence IDs.
    
    Args:
        fasta_file: Path to FASTA file
    
    Returns:
        Set of sequence IDs
    """
    ids = set()
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                # Extract ID (first word after >)
                seq_id = line[1:].strip().split()[0]
                ids.add(seq_id)
    return ids


def parse_interproscan_tsv(input_file: Path, fasta_ids: set = None, target_pfams: set = None) -> tuple:
    """
    Parse InterProScan TSV output and extract Pfam domain information.
    
    Args:
        input_file: Path to InterProScan TSV file
        fasta_ids: Set of IDs to include (from FASTA). If None, include all.
        target_pfams: Set of Pfam accessions to extract. If None, extract all Pfam domains.
    
    Returns:
        Tuple of (protein_domains dict, set of unique domain accessions)
    """
    protein_domains = defaultdict(lambda: {"size": 0, "domains": []})
    unique_domains = set()
    
    with open(input_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 8:
                protein_id = parts[0]
                protein_size = parts[2]
                database = parts[3]
                domain_acc = parts[4]
                start = parts[6]
                end = parts[7]
                
                # Skip if not in FASTA file (when fasta_ids provided)
                if fasta_ids is not None and protein_id not in fasta_ids:
                    continue
                
                protein_domains[protein_id]["size"] = protein_size
                
                # Check if it's a Pfam domain we're interested in
                if database == "Pfam":
                    if target_pfams is None or domain_acc in target_pfams:
                        protein_domains[protein_id]["domains"].append({
                            "acc": domain_acc,
                            "start": int(start),
                            "end": int(end)
                        })
                        unique_domains.add(domain_acc)
    
    return protein_domains, unique_domains


def generate_domain_config(unique_domains: set) -> dict:
    """
    Generate configuration (shape and color) for each unique domain.
    
    Args:
        unique_domains: Set of unique Pfam accessions
    
    Returns:
        Dictionary mapping Pfam accession to shape and color
    """
    domain_config = {}
    sorted_domains = sorted(unique_domains)
    
    for i, domain_acc in enumerate(sorted_domains):
        shape = SHAPES[i % len(SHAPES)]
        color = COLORS[i % len(COLORS)]
        domain_config[domain_acc] = {"shape": shape, "color": color, "name": domain_acc}
    
    return domain_config


def format_domain_definition(domain: dict, config: dict) -> str:
    """
    Format a domain into iTOL domain definition string.
    
    Args:
        domain: Dictionary with acc, start, end
        config: Domain configuration with shape and color
    
    Returns:
        Formatted string: SHAPE|start|end|color|PF#
    """
    return f"{config['shape']}|{domain['start']}|{domain['end']}|{config['color']}|{domain['acc']}"


def write_domain_table(protein_domains: dict, output_file: Path, domain_config: dict, fasta_ids: set = None):
    """
    Write domain summary table to TSV file.
    
    Args:
        protein_domains: Dictionary of protein domain information
        output_file: Path to output TSV file
        domain_config: Configuration for domain shapes and colors
        fasta_ids: Set of IDs from FASTA file (determines order and inclusion)
    """
    # Determine which IDs to output and in what order
    if fasta_ids is not None:
        # Use FASTA order (keep only those with domain info)
        output_ids = [pid for pid in sorted(fasta_ids) if pid in protein_domains]
    else:
        output_ids = sorted(protein_domains.keys())
    
    with open(output_file, "w") as out:
        out.write("Tree node ID\tTree node label\tProtein size\tDomain definitions\n")
        
        for pid in output_ids:
            info = protein_domains[pid]
            
            # Create label by replacing underscores with spaces
            label = pid.replace("_", " ")
            
            # Sort domains by start position (reversed - largest first)
            domains_sorted = sorted(info["domains"], key=lambda x: x["start"], reverse=True)
            
            # Format domain definitions
            domain_defs = []
            for d in domains_sorted:
                if d["acc"] in domain_config:
                    domain_defs.append(format_domain_definition(d, domain_config[d["acc"]]))
                else:
                    # Default formatting for unknown domains
                    default_config = {"shape": "RE", "color": "#808080"}
                    domain_defs.append(format_domain_definition(d, default_config))
            
            domain_str = "\t".join(domain_defs) if domain_defs else ""
            out.write(f"{pid}\t{label}\t{info['size']}\t{domain_str}\n")


def write_legend(domain_config: dict, output_file: Path):
    """
    Write a legend file showing domain colors and shapes.
    
    Args:
        domain_config: Dictionary mapping domains to shapes/colors
        output_file: Path to output legend file
    """
    legend_file = output_file.parent / f"{output_file.stem}_legend.tsv"
    with open(legend_file, "w") as out:
        out.write("Domain\tShape\tColor\n")
        for domain_acc in sorted(domain_config.keys()):
            config = domain_config[domain_acc]
            out.write(f"{domain_acc}\t{config['shape']}\t{config['color']}\n")
    print(f"Legend written to: {legend_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Create domain summary table from InterProScan TSV output"
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input InterProScan TSV file"
    )
    parser.add_argument(
        "-f", "--fasta",
        type=Path,
        required=True,
        help="FASTA file to filter proteins (only include proteins present in this file)"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output TSV file (default: INPUT_pfam_domains_summary.tsv)"
    )
    parser.add_argument(
        "-p", "--pfam",
        nargs="+",
        help="Specific Pfam accessions to extract (e.g., PF01134 PF21680). If not specified, extracts all Pfam domains."
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Optional config file for domain shapes and colors (TSV: pfam_acc, shape, color, name). "
             "If not provided, colors and shapes are auto-assigned."
    )
    parser.add_argument(
        "--legend",
        action="store_true",
        help="Generate a legend file showing domain colors and shapes"
    )
    
    args = parser.parse_args()
    
    # Check input files exist
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if not args.fasta.exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)
    
    # Set output file
    if args.output:
        output_file = args.output
    else:
        output_file = args.input.parent / f"{args.input.stem}_pfam_domains_summary.tsv"
    
    # Parse FASTA file to get IDs
    print(f"Reading FASTA file: {args.fasta}...")
    fasta_ids = parse_fasta_ids(args.fasta)
    print(f"Found {len(fasta_ids)} sequences in FASTA file")
    
    # Set target Pfam accessions
    target_pfams = set(args.pfam) if args.pfam else None
    
    # Parse domain file
    print(f"Parsing {args.input}...")
    protein_domains, unique_domains = parse_interproscan_tsv(args.input, fasta_ids, target_pfams)
    
    print(f"Found {len(unique_domains)} unique Pfam domains: {', '.join(sorted(unique_domains))}")
    
    # Generate or load domain configuration
    if args.config and args.config.exists():
        print(f"Loading domain config from: {args.config}")
        domain_config = {}
        with open(args.config) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    pfam_acc, shape, color = parts[0], parts[1], parts[2]
                    name = parts[3] if len(parts) > 3 else ""
                    domain_config[pfam_acc] = {"shape": shape, "color": color, "name": name}
        # Add auto-generated config for any domains not in config file
        for domain in unique_domains:
            if domain not in domain_config:
                auto_config = generate_domain_config({domain})
                domain_config.update(auto_config)
    else:
        print("Auto-generating domain colors and shapes...")
        domain_config = generate_domain_config(unique_domains)
    
    # Write output
    print(f"Writing {output_file}...")
    write_domain_table(protein_domains, output_file, domain_config, fasta_ids)
    
    # Write legend if requested
    if args.legend:
        write_legend(domain_config, output_file)
    
    # Summary
    proteins_in_output = sum(1 for pid in fasta_ids if pid in protein_domains)
    proteins_with_domains = sum(1 for pid in fasta_ids if pid in protein_domains and protein_domains[pid]["domains"])
    print(f"\nSummary:")
    print(f"  Proteins in FASTA: {len(fasta_ids)}")
    print(f"  Proteins with domain data: {proteins_in_output}")
    print(f"  Proteins with Pfam domains: {proteins_with_domains}")
    print("Done!")


if __name__ == "__main__":
    main()
