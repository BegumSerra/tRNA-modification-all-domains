import csv
import sys

# File paths from command line arguments
reference_tsv_path = sys.argv[1]  # TSV file containing the Pfam IDs
domain_results_tsv_path = sys.argv[2]   # TSV domain file to filter
output_tsv_path = sys.argv[3]  # Output FASTA file with matched entries

def extract_pfam_ids(ref_path):
    pfam_ids = set()
    with open(ref_path, 'r') as ref_file:
        ref_reader = csv.reader(ref_file, delimiter='\t')
        for row in ref_reader:
            if len(row) >= 5 and "Pfam" in row[3]:
                pfam_ids.add(row[4])
                print(pfam_ids)
    return pfam_ids

def create_domain_dict(domain_results_tsv_path):
    domain_dict = {}
    with open(domain_results_tsv_path, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            if len(row) >= 5 and "Pfam" in row[3]:
                key = row[0]
                domain_id = row[4]
                if key in domain_dict:
                    domain_dict[key].add(domain_id)
                else:
                    domain_dict[key] = {domain_id}
    print(domain_dict)
    return domain_dict

def filter_and_write_results(domain_results_tsv_path, domain_dict, pfam_ids, output_tsv_path):
    with open(domain_results_tsv_path, 'r') as infile, open(output_tsv_path, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        for row in reader:
            if len(row) >= 5 and "Pfam" in row[3]:
                key = row[0]
                domain_id = row[4]
                if key in domain_dict:
                    if domain_dict[key] == pfam_ids:
                        writer.writerow(row)


if __name__ == "__main__":
    # Extract Pfam IDs from the reference TSV file
    pfam_ids = extract_pfam_ids(reference_tsv_path)

    # Create a dictionary of domain results
    domain_dict = create_domain_dict(domain_results_tsv_path)

    # Filter the domain results dictionary and write matches to a new file
    filter_and_write_results(domain_results_tsv_path, domain_dict, pfam_ids, output_tsv_path)
