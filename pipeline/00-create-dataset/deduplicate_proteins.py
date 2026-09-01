"""
Deduplicate protein sequences by keeping only the first occurrence of each protein ID.
Updates both the protein FASTA file and the taxon_map.tsv file.

Exaample of duplicates found in the original FASTA file:
>WP_157042740.1 hypothetical protein [Chloroflexus aggregans]
>WP_157042740.1 hypothetical protein [Roseiflexus castenholzii]
"""

import argparse
from pathlib import Path


def deduplicate_fasta(input_fasta, output_fasta):
    """
    Remove duplicate protein IDs from FASTA file, keeping only first occurrence.
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output deduplicated FASTA file
    
    Returns:
        tuple: (set of unique IDs kept, set of duplicate IDs removed)
    """
    seen_ids = set()
    duplicate_ids = set()
    kept_sequences = []
    
    current_header = None
    current_sequence = []
    
    with open(input_fasta, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            
            if line.startswith('>'):
                # Save previous sequence if it was unique
                if current_header is not None:
                    protein_id = current_header.split()[0].lstrip('>')
                    if protein_id not in seen_ids:
                        seen_ids.add(protein_id)
                        kept_sequences.append((current_header, ''.join(current_sequence)))
                    else:
                        duplicate_ids.add(protein_id)
                
                # Start new sequence
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Don't forget the last sequence
        if current_header is not None:
            protein_id = current_header.split()[0].lstrip('>')
            if protein_id not in seen_ids:
                seen_ids.add(protein_id)
                kept_sequences.append((current_header, ''.join(current_sequence)))
            else:
                duplicate_ids.add(protein_id)
    
    # Write deduplicated sequences
    with open(output_fasta, 'w') as f:
        for header, sequence in kept_sequences:
            f.write(f"{header}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")
    
    return seen_ids, duplicate_ids


def update_taxon_map(input_taxon_map, output_taxon_map, kept_ids):
    """
    Update taxon_map to only include protein IDs that were kept.
    
    Args:
        input_taxon_map: Path to input taxon_map.tsv
        output_taxon_map: Path to output deduplicated taxon_map.tsv
        kept_ids: Set of protein IDs to keep
    """
    kept_count = 0
    removed_count = 0
    
    with open(input_taxon_map, 'r') as f_in, open(output_taxon_map, 'w') as f_out:
        for line in f_in:
            protein_id = line.split('\t')[0]
            if protein_id in kept_ids:
                f_out.write(line)
                kept_count += 1
            else:
                removed_count += 1
    
    return kept_count, removed_count


def main():
    parser = argparse.ArgumentParser(
        description='Deduplicate protein sequences and update taxon_map'
    )
    parser.add_argument(
        '--input-fasta',
        required=True,
        help='Input protein FASTA file'
    )
    parser.add_argument(
        '--output-fasta',
        required=True,
        help='Output deduplicated FASTA file'
    )
    parser.add_argument(
        '--input-taxon-map',
        required=True,
        help='Input taxon_map.tsv file'
    )
    parser.add_argument(
        '--output-taxon-map',
        required=True,
        help='Output deduplicated taxon_map.tsv file'
    )
    
    args = parser.parse_args()
    
    print(f"Deduplicating {args.input_fasta}...")
    kept_ids, duplicate_ids = deduplicate_fasta(args.input_fasta, args.output_fasta)
    
    print(f"  Kept: {len(kept_ids)} unique protein IDs")
    print(f"  Removed: {len(duplicate_ids)} duplicate protein IDs")
    
    if duplicate_ids:
        print(f"  Example duplicates: {', '.join(list(duplicate_ids)[:5])}")
    
    print(f"\nUpdating {args.input_taxon_map}...")
    kept_count, removed_count = update_taxon_map(
        args.input_taxon_map,
        args.output_taxon_map,
        kept_ids
    )
    
    print(f"  Kept: {kept_count} taxon_map entries")
    print(f"  Removed: {removed_count} taxon_map entries")
    
    print("\nDeduplication complete!")


if __name__ == '__main__':
    main()
