#!/usr/bin/env python3
"""
Download protein sequences from NCBI for a list of assembly accessions
Works for both bacteria and archaea
"""
import subprocess
import sys
from pathlib import Path
import argparse


def download_proteins_batch(accessions, output_dir, batch_size=50):
    """Download proteins in batches using NCBI datasets tool"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    total = len(accessions)
    num_batches = (total + batch_size - 1) // batch_size
    
    for i in range(0, total, batch_size):
        batch_num = i // batch_size + 1
        batch_accessions = accessions[i:i + batch_size]
        batch_dir = output_dir / f"batch_{batch_num}"
        
        print(f"\n{'='*60}")
        print(f"Batch {batch_num}/{num_batches}: Downloading {len(batch_accessions)} assemblies")
        print(f"{'='*60}")
        
        # Create batch directory
        batch_dir.mkdir(exist_ok=True)
        
        # Download using datasets command
        cmd = [
            'datasets', 'download', 'genome', 'accession',
            *batch_accessions,
            '--include', 'protein',
            '--filename', str(batch_dir / 'data.zip')
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"✓ Downloaded batch {batch_num}")
            
            # Unzip
            unzip_cmd = ['unzip', '-q', '-o', str(batch_dir / 'data.zip'), '-d', str(batch_dir)]
            subprocess.run(unzip_cmd, check=True)
            print(f"✓ Extracted batch {batch_num}")
            
        except subprocess.CalledProcessError as e:
            print(f"✗ Error downloading batch {batch_num}: {e}")
            print(f"  stderr: {e.stderr}")
            continue
    
    print(f"\n{'='*60}")
    print(f"Download complete!")
    print(f"{'='*60}\n")


def combine_proteins(output_dir, combined_file):
    """Combine all downloaded protein files into one"""
    output_dir = Path(output_dir)
    combined_file = Path(combined_file)
    
    print(f"\nCombining protein files into {combined_file}...")
    
    # Find all protein.faa files
    protein_files = list(output_dir.glob("batch_*/ncbi_dataset/data/GC*/protein.faa"))
    
    if not protein_files:
        print("Warning: No protein.faa files found!")
        return
    
    print(f"Found {len(protein_files)} protein files")
    
    # Combine files
    combined_file.parent.mkdir(parents=True, exist_ok=True)
    with open(combined_file, 'w') as outf:
        for prot_file in sorted(protein_files):
            with open(prot_file, 'r') as inf:
                outf.write(inf.read())
    
    print(f"✓ Combined {len(protein_files)} files into {combined_file}")
    
    # Count sequences
    with open(combined_file, 'r') as f:
        seq_count = sum(1 for line in f if line.startswith('>'))
    print(f"✓ Total protein sequences: {seq_count:,}")


def main():
    parser = argparse.ArgumentParser(description='Download protein sequences from NCBI')
    parser.add_argument('accession_list', help='File with assembly accessions (one per line)')
    parser.add_argument('output_dir', help='Output directory for downloads')
    parser.add_argument('combined_output', help='Output file for combined proteins')
    parser.add_argument('--batch-size', type=int, default=50, help='Batch size for downloads (default: 50)')
    
    args = parser.parse_args()
    
    # Read accession list
    print(f"Reading accessions from {args.accession_list}...")
    with open(args.accession_list, 'r') as f:
        accessions = [line.strip() for line in f if line.strip() and not line.startswith('#') and line.strip() != 'assembly_accession']
    
    print(f"Found {len(accessions)} accessions to download\n")
    
    # Download proteins
    download_proteins_batch(accessions, args.output_dir, args.batch_size)
    
    # Combine into single file
    combine_proteins(args.output_dir, args.combined_output)
    
    print("\n✓ All done!")


if __name__ == '__main__':
    main()
