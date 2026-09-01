import sys
import os
from pathlib import Path

def extract_euk_seqs(input_fasta):
    """Extracts eukaryotic sequences (headers starting with 'EP' and numbers) from a fasta file."""
    euk_seqs = []
    with open(input_fasta) as f:
        write_seq = False
        seq_lines = []
        header = ''
        for line in f:
            if line.startswith('>'):
                if seq_lines and write_seq:
                    euk_seqs.append((header, ''.join(seq_lines)))
                header = line.strip()
                # Check if header starts with >EP and numbers
                if header[1:].startswith('EP') and header[3:8].isdigit():
                    write_seq = True
                else:
                    write_seq = False
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        # Add last sequence if needed
        if seq_lines and write_seq:
            euk_seqs.append((header, ''.join(seq_lines)))
    return euk_seqs

def main():
    if len(sys.argv) < 3:
        print("Usage: python extract-euk-fa.py <input_fasta> <output_fasta>")
        sys.exit(1)
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    euk_seqs = extract_euk_seqs(input_fasta)
    with open(output_fasta, 'w') as out:
        for header, seq in euk_seqs:
            out.write(f"{header}\n")
            # Write sequence in lines of 60 chars
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')
    print(f"Extracted {len(euk_seqs)} eukaryotic sequences to {output_fasta}")

if __name__ == "__main__":
    main()
