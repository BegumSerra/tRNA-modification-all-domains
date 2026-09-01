import csv
import sys

# Assuming command line arguments for file paths
domain_tsv_path = sys.argv[1]  # Domain TSV file to extract IDs from
blast_result_path = sys.argv[2]  # BLAST result file to filter
output_file_path = sys.argv[3]  # Output file for filtered rows


def normalize_id(raw_id):
    """Normalize an ID by removing prefixes/suffixes and returning the base ID."""
    # This example splits by '|' and takes the middle part, adjust as necessary
    parts = raw_id.split('|')
    return parts[-2] if len(parts) > 1 else parts[0]


def extract_ids(domain_tsv_path):
    """Extract and normalize IDs from the first column of the domain TSV file."""
    ids = set()
    with open(domain_tsv_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row:  # Ensure the row is not empty
                normalized_id = normalize_id(row[0])
                ids.add(normalized_id)
    return ids


def filter_blast_results(blast_result_path, ids, output_file_path):
    """Filter the BLAST result file using normalized IDs and write matching rows to the output file."""
    with open(blast_result_path, 'r') as infile, open(output_file_path, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for row in reader:
            if len(row) > 1:
                normalized_blast_id = normalize_id(row[1])
                if normalized_blast_id in ids:
                    writer.writerow(row)


if __name__ == "__main__":
    # Extract and normalize IDs from the domain file
    ids = extract_ids(domain_tsv_path)

    # Use the normalized IDs to filter the BLAST result file
    filter_blast_results(blast_result_path, ids, output_file_path)
