"""
Merge filter IDs and exclude those already in KEGG filter list
"""
import sys


def load_ids(file_path):
    """Load IDs from a file into a set"""
    try:
        with open(file_path, 'r') as f:
            return set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        return set()


def main(wo_domain_file, w_extra_domain_file, kegg_filter_file, output_file):
    # Load IDs from domain filter files
    wo_domain_ids = load_ids(wo_domain_file)
    w_extra_domain_ids = load_ids(w_extra_domain_file)
    
    # Combine both sets
    all_domain_filter_ids = wo_domain_ids | w_extra_domain_ids
    
    # Load KEGG filter IDs
    kegg_filter_ids = load_ids(kegg_filter_file)
    
    # Remove IDs that are already in KEGG filter
    final_ids = all_domain_filter_ids - kegg_filter_ids
    
    print(f"Domain filter IDs (wo_domain): {len(wo_domain_ids)}")
    print(f"Domain filter IDs (w_extra_domain): {len(w_extra_domain_ids)}")
    print(f"Combined domain filter IDs: {len(all_domain_filter_ids)}")
    print(f"KEGG filter IDs: {len(kegg_filter_ids)}")
    print(f"IDs already in KEGG filter (removed from domain filter): {len(kegg_filter_ids & all_domain_filter_ids)}")
    print(f"Final IDs to filter (after excluding KEGG duplicates): {len(final_ids)}")
    
    # Save final IDs
    with open(output_file, 'w') as f:
        for id in sorted(final_ids):
            f.write(id + '\n')
    
    print(f"Final filter IDs saved to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python merge-and-exclude-kegg.py <wo_domain_file> <w_extra_domain_file> <kegg_filter_file> <output_file>")
        sys.exit(1)
    
    wo_domain_file = sys.argv[1]
    w_extra_domain_file = sys.argv[2]
    kegg_filter_file = sys.argv[3]
    output_file = sys.argv[4]
    
    main(wo_domain_file, w_extra_domain_file, kegg_filter_file, output_file)
