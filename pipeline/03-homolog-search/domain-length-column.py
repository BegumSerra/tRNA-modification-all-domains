import csv
import sys

# File paths from command line arguments
domain_tsv_path = sys.argv[1]  # The path to the TSV file containing domain information
output_csv_path = sys.argv[2]  # The path for the output file with the new "Domain Length" column

def process_file(domain_tsv_path, output_csv_path):
    with open(domain_tsv_path, 'r') as infile, open(output_csv_path, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Optional: Write header row to the output file, adjust according to your specific needs
        # writer.writerow(['Your', 'Header', 'Columns', 'Including', 'New', 'Domain Length'])

        for row in reader:
            if len(row) > 7 and "Pfam" in row[3]:  # Check if the row meets criteria and has enough columns
                start_position = int(row[6])
                stop_position = int(row[7])
                domain_length = stop_position - start_position + 1
                row.insert(6, domain_length)  # Insert domain length before the 'Start' column
            writer.writerow(row)

if __name__ == "__main__":
    process_file(domain_tsv_path, output_csv_path)
