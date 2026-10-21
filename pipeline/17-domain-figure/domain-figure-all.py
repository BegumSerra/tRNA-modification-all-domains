# python
import argparse
import os
import matplotlib
matplotlib.use("Agg")  # non-interactive backend; avoids Qt/display dependency issues
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import yaml


def read_config(config_file):
    """Read config file with protein definitions

    YAML format:
    proteins:
      - name: ArcS
        pfam: /path/to/ArcS.tsv
        fasta: /path/to/ArcS.fasta
        protein_id: WP_248517187.1
      - name: MnmA
        pfam: /path/to/MnmA.tsv
        fasta: /path/to/MnmA.fasta
        protein_id: WP_001234567.1

    Or simple TSV format (tab-separated):
    # name	pfam	fasta	protein_id
    ArcS	/path/to/ArcS.tsv	/path/to/ArcS.fasta	WP_248517187.1
    MnmA	/path/to/MnmA.tsv	/path/to/MnmA.fasta	WP_001234567.1
    """
    if config_file.endswith(('.yaml', '.yml')):
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config.get('proteins', [])
    else:
        # TSV format
        proteins = []
        with open(config_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    proteins.append({
                        'name': parts[0],
                        'pfam': parts[1],
                        'fasta': parts[2],
                        'protein_id': parts[3]
                    })
        return proteins


def read_pfam_domains(pfam_tsv):
    """Read Pfam TSV and extract relevant columns"""
    try:
        df = pd.read_csv(pfam_tsv, sep="\t", header=None, dtype=str, engine='python')
    except Exception as e:
        print(f"Warning: Error reading with python engine: {e}")
        rows = []
        with open(pfam_tsv, 'r') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                rows.append(parts)
        max_cols = max(len(row) for row in rows) if rows else 0
        for row in rows:
            while len(row) < max_cols:
                row.append('')
        df = pd.DataFrame(rows)

    if df.shape[1] >= 8:
        pfam_mask = df[3].astype(str).str.contains('Pfam', case=False, na=False)
        df_pfam = df[pfam_mask].copy()

        if len(df_pfam) == 0:
            print(f"Warning: No Pfam entries found in {pfam_tsv}")
            return pd.DataFrame(columns=['protein_id', 'pfam_id', 'start', 'end'])

        out = pd.DataFrame({
            'protein_id': df_pfam[0].astype(str),
            'pfam_id': df_pfam[4].astype(str),
            'start': pd.to_numeric(df_pfam[6], errors='coerce'),
            'end': pd.to_numeric(df_pfam[7], errors='coerce')
        })
    else:
        raise ValueError(f"Pfam TSV has {df.shape[1]} columns, expected at least 8")

    out = out.dropna(subset=['start', 'end']).copy()
    out['start'] = out['start'].astype('Int64')
    out['end'] = out['end'].astype('Int64')
    return out.reset_index(drop=True)


def read_blast_lengths(blast_tsv):
    df = pd.read_csv(blast_tsv, sep="\t", header=None, dtype=str)
    protein_ids = df[1].astype(str).str.replace(r'^[a-z]+\|', '', regex=True).str.replace(r'\|$', '', regex=True)
    out = pd.DataFrame({
        "protein_id": protein_ids,
        "length": pd.to_numeric(df[3], errors="coerce").astype("Int64"),
        "species_name": ""
    })
    out = out.dropna(subset=["length"]).copy()
    out["length"] = out["length"].astype(int)
    return out.reset_index(drop=True)


def read_fasta_lengths(fasta_file):
    """Read protein lengths from FASTA file"""
    proteins = []
    current_id = None
    current_species = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if current_id is not None:
                    seq_length = len(''.join(current_seq))
                    proteins.append({'protein_id': current_id, 'length': seq_length, 'species_name': current_species or ''})

                header = line[1:].strip()
                species_name = ""
                if '[' in header and ']' in header:
                    start = header.rfind('[')
                    end = header.rfind(']')
                    if start < end:
                        species_name = header[start+1:end].strip()

                protein_id = header.split()[0]
                protein_id = protein_id.replace('ref|', '').replace('sp|', '').replace('gb|', '').replace('dbj|', '').rstrip('|')
                current_id = protein_id
                current_species = species_name
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            seq_length = len(''.join(current_seq))
            proteins.append({'protein_id': current_id, 'length': seq_length, 'species_name': current_species or ''})

    if not proteins:
        raise ValueError(f"No sequences found in FASTA file: {fasta_file}")

    return pd.DataFrame(proteins).reset_index(drop=True)


def sanitize_name(name):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(name))



def plot_multiple_protein_domains(proteins_data, domains_df, outpath, dpi=200):
    """Plot multiple proteins stacked vertically - each on separate line"""
    num_proteins = len(proteins_data)
    if num_proteins == 0:
        return

    # Figure setup
    fig_height = num_proteins * 0.3 + 0.5
    fig, axes = plt.subplots(num_proteins, 1, figsize=(18, fig_height), dpi=dpi)

    if num_proteins == 1:
        axes = [axes]

    # Get max length for scaling
    max_length = max([p['length'] for p in proteins_data])

    # Create GLOBAL color map for ALL Pfam IDs across all proteins
    # This ensures same Pfam ID has same color in all proteins
    all_protein_ids = [p['protein_id'] for p in proteins_data]
    all_domains = domains_df[domains_df['protein_id'].isin(all_protein_ids)]
    all_pfams = all_domains['pfam_id'].unique().tolist()
    cmap = plt.get_cmap('Pastel2')
    global_colors = {p: cmap(i / max(len(all_pfams), 1)) for i, p in enumerate(all_pfams)}

    for idx, (ax, pinfo) in enumerate(zip(axes, proteins_data)):
        pid = pinfo['protein_id']
        length = int(pinfo['length'])
        pname = pinfo.get('protein_name', '')
        species = pinfo.get('species_name', '')

        # Get domains for THIS protein only
        doms = domains_df[domains_df['protein_id'] == pid].copy()
        doms = doms.drop_duplicates(subset=['pfam_id', 'start', 'end'])
        doms = doms.sort_values('start')

        print(f"    {pid}: {len(doms)} domains, {length} aa")

        # Setup axis
        ax.set_xlim(0, max_length + 50)
        ax.set_ylim(0, 1)
        ax.axis('off')
        y = 0.5

        # Draw protein name and species on left (further left with bigger margin)
        ax.text(-0.22, 0.7, pname, transform=ax.transAxes, fontsize=12, weight='bold', ha='left')
        ax.text(-0.22, 0.1, species, transform=ax.transAxes, fontsize=9, style='italic', ha='left', color='#333')

        # Draw baseline
        ax.hlines(y, 0, length, color='black', linewidth=3, zorder=1)
        ax.text(length + 5, y, f"{length} aa", fontsize=9, va='center')

        # Draw domains using global color map
        rect_height = 1
        for _, r in doms.iterrows():
            s, e = int(r['start']), int(r['end'])
            width = e - s
            # Use round boxstyle with mutation_scale for squircle-like corners
            box = FancyBboxPatch((s, y - rect_height/2), width, rect_height,
                                 boxstyle="round,pad=0.01",
                                 facecolor=global_colors.get(r['pfam_id'], (0.7,0.7,0.7)),
                                 edgecolor='black', linewidth=0, zorder=2,
                                 mutation_scale=15)
            ax.add_patch(box)
            ax.text(s + width/2, y, r['pfam_id'], ha='center', va='center', fontsize=12, weight='bold', zorder=3, fontname='Arial')

    plt.tight_layout()
    fig.subplots_adjust(left=0.25)  # Bigger left margin
    fig.savefig(outpath, format='svg', bbox_inches='tight')
    plt.close(fig)


def main(pfam_files=None, protein_files=None, outdir=None, protein_ids=None, file_type='auto', config_file=None):
    os.makedirs(outdir, exist_ok=True)

    all_pfam_data = pd.DataFrame()
    all_proteins_data = []

    # If config file is provided, use it
    if config_file:
        print(f"Reading config from: {config_file}")
        config_proteins = read_config(config_file)

        for item in config_proteins:
            pfam_tsv = item['pfam']
            protein_file = item.get('fasta') or item.get('blast')
            protein_name = item.get('name', os.path.splitext(os.path.basename(pfam_tsv))[0])
            pid = item['protein_id']

            print(f"\nProcessing: {protein_name} - {pid}")
            print(f"  Pfam: {pfam_tsv}")
            print(f"  Proteins: {protein_file}")

            # Read pfam data
            pfam = read_pfam_domains(pfam_tsv)
            print(f"  Loaded {len(pfam)} domain records")

            # Auto-detect file type
            ftype = file_type
            if ftype == 'auto':
                if protein_file.endswith(('.fasta', '.fa', '.faa')):
                    ftype = 'fasta'
                else:
                    ftype = 'blast'

            # Read protein lengths
            if ftype == 'fasta':
                proteins = read_fasta_lengths(protein_file)
            else:
                proteins = read_blast_lengths(protein_file)

            # Combine pfam data
            all_pfam_data = pd.concat([all_pfam_data, pfam], ignore_index=True)

            # Get protein data
            row = proteins[proteins['protein_id'] == pid]
            if len(row) > 0:
                all_proteins_data.append({
                    'protein_id': pid,
                    'length': int(row.iloc[0]['length']),
                    'protein_name': protein_name,
                    'species_name': str(row.iloc[0].get('species_name', ''))
                })
                print(f"  ✓ Found protein: {pid} ({row.iloc[0]['length']} aa)")
            else:
                print(f"  ✗ Protein {pid} not found in {protein_file}")

    else:
        # Original logic for command-line arguments
        if isinstance(pfam_files, str):
            pfam_files = [pfam_files]
        if isinstance(protein_files, str):
            protein_files = [protein_files]

        # Process each file pair
        for i, (pfam_tsv, protein_file) in enumerate(zip(
            pfam_files * len(protein_files) if len(pfam_files) == 1 else pfam_files,
            protein_files * len(pfam_files) if len(protein_files) == 1 else protein_files
        )):
            protein_name = os.path.splitext(os.path.basename(pfam_tsv))[0]

            print(f"\nReading Pfam domains from: {pfam_tsv}")
            pfam = read_pfam_domains(pfam_tsv)
            print(f"  Loaded {len(pfam)} domain records for {pfam['protein_id'].nunique()} proteins")

            ftype = file_type
            if ftype == 'auto':
                if protein_file.endswith(('.fasta', '.fa', '.faa')):
                    ftype = 'fasta'
                else:
                    ftype = 'blast'

            print(f"Reading protein file ({ftype}): {protein_file}")
            if ftype == 'fasta':
                proteins = read_fasta_lengths(protein_file)
            else:
                proteins = read_blast_lengths(protein_file)
            print(f"  Loaded {len(proteins)} proteins")

            all_pfam_data = pd.concat([all_pfam_data, pfam], ignore_index=True)

            if protein_ids:
                pids = protein_ids if i == 0 else []
                if len(pfam_files) > 1 and len(protein_ids) > i:
                    pids = [protein_ids[i]] if isinstance(protein_ids[i], str) else protein_ids[i]
            else:
                print("\nAvailable proteins:")
                for idx, pid in enumerate(proteins['protein_id'].unique()[:20], 1):
                    print(f"  {idx}. {pid}")
                choice = input("Enter protein ID(s) comma-separated: ").strip()
                pids = [p.strip() for p in choice.split(',')]

            for pid in pids:
                row = proteins[proteins['protein_id'] == pid]
                if len(row) > 0:
                    all_proteins_data.append({
                        'protein_id': pid,
                        'length': int(row.iloc[0]['length']),
                        'protein_name': protein_name,
                        'species_name': str(row.iloc[0].get('species_name', ''))
                    })

    if len(all_proteins_data) == 0:
        print("No matching proteins found!")
        return

    print(f"\nGenerating figure for {len(all_proteins_data)} protein(s)...")

    if len(all_proteins_data) > 1:
        outfile = os.path.join(outdir, "combined_stacked_TadA_ADAT2.svg")
    else:
        outfile = os.path.join(outdir, f"{sanitize_name(all_proteins_data[0]['protein_id'])}.svg")

    plot_multiple_protein_domains(all_proteins_data, all_pfam_data, outfile)

    print(f"\n✓ Saved: {outfile}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Pfam domain architecture")
    parser.add_argument("--config", help="Config file (YAML or TSV) with protein definitions")
    parser.add_argument("--pfam", nargs='+', help="Pfam TSV file(s)")
    parser.add_argument("--proteins", nargs='+', help="FASTA or BLAST file(s)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--protein-id", default=None, nargs='+', help="Protein ID(s)")
    parser.add_argument("--file-type", default='auto', choices=['auto', 'fasta', 'blast'])
    args = parser.parse_args()

    # Check that either config or pfam+proteins is provided
    if not args.config and not (args.pfam and args.proteins):
        parser.error("Either --config or both --pfam and --proteins are required")

    # Parse protein IDs
    protein_ids = None
    if args.protein_id:
        protein_ids = []
        for pid in args.protein_id:
            protein_ids.extend([p.strip() for p in pid.split(',')])

    main(
        pfam_files=args.pfam,
        protein_files=args.proteins,
        outdir=args.outdir,
        protein_ids=protein_ids,
        file_type=args.file_type,
        config_file=args.config
    )

