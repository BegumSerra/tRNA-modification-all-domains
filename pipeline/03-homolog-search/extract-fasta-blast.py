from Bio import SeqIO
import sys

# Specify your file paths here
blast_file_path = sys.argv[1]
proteome_file_path = sys.argv[2]
output_file_path = sys.argv[3]

def extract_ids_from_blast(blast_file_path):
    ids = []
    with open(blast_file_path, 'r') as blast_file:
        for line in blast_file:
            parts = line.strip().split('\t')  # Ensure to strip newline characters
            # Check if there are enough parts in the line
            if len(parts) > 1:
                id_full = parts[1]
                # Extract the inner part if it starts with 'gb|' or 'ref|', else use the whole ID
                if id_full.startswith('gb|') or id_full.startswith('ref|'):
                    id_part = id_full.split('|')[1]
                    ids.append(id_part)
                else:
                    # For IDs without 'gb|' or 'ref|', add the whole ID
                    ids.append(id_full)
    return ids


def extract_fasta_sequences(proteome_file_path, ids, output_file_path):
    sequences = []
    for record in SeqIO.parse(proteome_file_path, "fasta"):
        if any(id_ in record.id for id_ in ids):
            sequences.append(record)
    SeqIO.write(sequences, output_file_path, "fasta")
    print(f"Extracted {len(sequences)} sequences to {output_file_path}")

# Extract IDs from the BLAST file
ids = extract_ids_from_blast(blast_file_path)

# Extract matching FASTA sequences from the proteome file and save them
extract_fasta_sequences(proteome_file_path, ids, output_file_path)
