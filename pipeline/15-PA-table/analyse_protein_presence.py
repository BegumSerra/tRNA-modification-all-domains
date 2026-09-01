#!/usr/bin/env python3
import os
import sys
import pandas as pd
from collections import defaultdict
import glob
import re

def count_sequences_in_fasta(fasta_file):
    """Count the number of sequences in a fasta file."""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count

def extract_protein_id_from_header(header):
    """Extract protein ID (like WP_046437935.1) from fasta header."""
    # Look for patterns like WP_123456789.1, NP_123456789.1, etc.
    protein_id_match = re.search(r'([A-Z]{2}_\d+\.\d+)', header)
    if protein_id_match:
        return protein_id_match.group(1)
    return None

def extract_species_from_header(header):
    """Extract species name from fasta header with multiple strategies."""
    # Remove '>' and clean up
    header = header.strip().lstrip('>')
    
    # Strategy 1: Handle MAST cases - EP00944_MAST-04C_sp_MAST-4C-sp1_P004422
    # Extract and convert to format: MAST-04C sp MAST-4C-sp1
    mast_pattern = re.search(r'EP\d+_(MAST-\d+[A-Z])_sp_(MAST-\d+[A-Z]-sp\d+)', header)
    if mast_pattern:
        part1 = mast_pattern.group(1)  # MAST-04C
        part2 = mast_pattern.group(2)  # MAST-4C-sp1
        return f"{part1} sp {part2}"
    
    # Strategy 2: Handle EP cases - extract everything between EP####_ and _P######
    # This handles: EP00742_Colponemidia_sp_Colp-10_P027420 -> Colponemidia sp Colp-10
    # And also: EP00682_Trypanosoma_cruzi_P010917 -> Trypanosoma cruzi
    ep_full_pattern = re.search(r'EP\d+_(.+?)_P\d+', header)
    if ep_full_pattern:
        species_with_underscore = ep_full_pattern.group(1)
        return species_with_underscore.replace('_', ' ')
    
    # Strategy 3: Look for species name in square brackets [Genus species] - for ALL other cases
    bracket_match = re.search(r'\[([^\]]+)\]', header)
    if bracket_match:
        species_name = bracket_match.group(1).strip()
        # Handle cases like "Candidatus Genus species"
        parts = species_name.split()
        if len(parts) >= 2:
            # If it starts with "Candidatus", take the first two words after "Candidatus"
            if parts[0].lower() == 'candidatus' and len(parts) >= 3:
                return f"{parts[1]} {parts[2]}"  # Take genus and species after Candidatus
            else:
                return f"{parts[0]} {parts[1]}"  # Normal genus species
        return species_name
    
    # Strategy 4: Look for OS=Genus species pattern
    os_match = re.search(r'OS=([A-Za-z][A-Za-z_\-\d]*\s+[A-Za-z][A-Za-z_\-\d]*)', header)
    if os_match:
        return os_match.group(1)
    
    # Strategy 5: Only as last resort - look for patterns like "Genus species" at the beginning
    # But skip if it looks like an accession number
    parts = header.split()
    if len(parts) >= 2:
        # Skip if first part looks like an accession number (contains numbers and dots)
        if not re.search(r'[0-9\.]', parts[0]):
            # Check if first two parts look like genus species
            if re.match(r'^[A-Za-z][A-Za-z_\-\d]*$', parts[0]) and re.match(r'^[A-Za-z][A-Za-z_\-\d]*$', parts[1]):
                return f"{parts[0]} {parts[1]}"
    
    return ""

def extract_species_from_ncbi_taxonomy(ncbi_taxonomy):
    """Extract species name from NCBI taxonomy string."""
    if not ncbi_taxonomy or ncbi_taxonomy == 'nan':
        return None
    
    # NCBI taxonomy format: d__;p__;c__;o__;f__;g__;s__
    parts = ncbi_taxonomy.split(';')
    
    for part in parts:
        part = part.strip()
        if part.startswith('s__'):
            species = part.replace('s__', '').strip()
            if species and species != '':
                return species
    
    return None

def normalize_species_name(name):
    """Normalize species names for better matching."""
    if not name:
        return ""
    
    # Remove extra spaces and convert to lowercase for comparison
    name = re.sub(r'\s+', ' ', name.strip())
    
    # Handle common variations
    name = re.sub(r'_[A-Z]$', '', name)  # Remove suffix like _A, _B
    name = re.sub(r'_[A-Z]\s', ' ', name)  # Remove _A from middle
    
    return name

def handle_candidatus_special_case(species_name, tsv_species_list):
    """Handle special case for Candidatus species - try matching with and without 'Candidatus'."""
    if not species_name:
        return None
    
    species_lower = species_name.lower()
    
    # Create normalized TSV species lookup
    normalized_tsv_species = {normalize_species_name(s).lower(): s for s in tsv_species_list}
    
    # If the extracted species doesn't contain "candidatus", try adding it
    if 'candidatus' not in species_lower:
        # Try with "Candidatus" prefix
        candidatus_version = f"candidatus {species_name.lower()}"
        if candidatus_version in normalized_tsv_species:
            return normalized_tsv_species[candidatus_version]
    
    # If the extracted species contains "candidatus", try without it
    if 'candidatus' in species_lower:
        # Remove "candidatus" and try matching
        non_candidatus_version = species_lower.replace('candidatus', '').strip()
        if non_candidatus_version in normalized_tsv_species:
            return normalized_tsv_species[non_candidatus_version]
    
    # Try partial matching for both versions
    genus_species = species_name.lower().replace('candidatus', '').strip()
    parts = genus_species.split()
    
    if len(parts) >= 2:
        genus = parts[0]
        species_epithet = parts[1]
        
        # Try various combinations
        combinations_to_try = [
            f"candidatus {genus} {species_epithet}",
            f"{genus} {species_epithet}",
            species_lower,
            genus_species
        ]
        
        for combo in combinations_to_try:
            if combo in normalized_tsv_species:
                return normalized_tsv_species[combo]
        
        # Try partial matching
        for norm_tsv_name, orig_tsv_name in normalized_tsv_species.items():
            # Check if genus and species epithet are both present
            if genus in norm_tsv_name and species_epithet in norm_tsv_name:
                return orig_tsv_name
    
    return None

def load_protein_genome_mapping(mapping_file):
    """Load the protein genome mapping file for additional species matching."""
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        print(f"Loaded protein genome mapping with {len(mapping_df)} entries")
        
        # Create a mapping dictionary - we'll use multiple columns for matching
        genome_to_species = {}
        protein_id_to_species = {}
        
        # Check what columns are available
        print(f"Mapping file columns: {mapping_df.columns.tolist()}")
        
        # Create various mappings based on available columns
        for _, row in mapping_df.iterrows():
            # Try different column combinations for species identification
            potential_species = []
            
            # Check for organism or species columns
            for col in mapping_df.columns:
                col_lower = col.lower()
                if any(keyword in col_lower for keyword in ['organism', 'species', 'name', 'binomial']):
                    if pd.notna(row[col]) and str(row[col]).strip():
                        potential_species.append(str(row[col]).strip())
            
            # Also try to extract from NCBI taxonomy
            ncbi_taxonomy = str(row.get('ncbi_taxonomy', ''))
            if ncbi_taxonomy and ncbi_taxonomy != 'nan':
                species_from_tax = extract_species_from_ncbi_taxonomy(ncbi_taxonomy)
                if species_from_tax:
                    potential_species.append(species_from_tax)
            
            # Check for genome_id or accession
            genome_id = None
            for col in mapping_df.columns:
                col_lower = col.lower()
                if any(keyword in col_lower for keyword in ['genome', 'accession', 'id']):
                    if pd.notna(row[col]) and str(row[col]).strip():
                        genome_id = str(row[col]).strip()
                        break
            
            # Check for protein_id column AND any protein IDs in the row
            protein_ids = []
            for col in mapping_df.columns:
                col_lower = col.lower()
                if 'protein' in col_lower and ('id' in col_lower or 'accession' in col_lower):
                    if pd.notna(row[col]) and str(row[col]).strip():
                        protein_ids.append(str(row[col]).strip())
                
                # Also check if any cell contains protein ID patterns
                cell_value = str(row[col])
                if pd.notna(row[col]) and cell_value != 'nan':
                    # Look for protein ID patterns in any column
                    protein_id_matches = re.findall(r'([A-Z]{2}_\d+\.\d+)', cell_value)
                    protein_ids.extend(protein_id_matches)
            
            # Store all potential mappings
            for species in potential_species:
                if species and len(species.split()) >= 2:  # Ensure it's a proper binomial name
                    normalized = normalize_species_name(species).lower()
                    genome_to_species[normalized] = species
                    
                    # Also store with genome_id if available
                    if genome_id:
                        genome_to_species[genome_id.lower()] = species
                    
                    # Store with all protein_ids if available
                    for protein_id in protein_ids:
                        if protein_id:
                            protein_id_to_species[protein_id] = species
        
        print(f"Created {len(genome_to_species)} species mappings from protein genome file")
        print(f"Created {len(protein_id_to_species)} protein ID mappings")
        
        # Store the full dataframe for direct protein ID lookups
        return genome_to_species, protein_id_to_species, mapping_df
        
    except Exception as e:
        print(f"Error loading protein genome mapping: {e}")
        return {}, {}, pd.DataFrame()

def lookup_protein_id_in_mapping(protein_id, mapping_df, tsv_species_list):
    """Look up a protein ID directly in the mapping dataframe and return matched species."""
    if protein_id and not mapping_df.empty:
        # Search for the protein ID in any column of the mapping file
        for _, row in mapping_df.iterrows():
            # Check each column for the protein ID
            for col in mapping_df.columns:
                cell_value = str(row[col])
                if protein_id in cell_value:
                    # Found the protein ID, now extract species from this row
                    # Try NCBI taxonomy first
                    ncbi_taxonomy = str(row.get('ncbi_taxonomy', ''))
                    if ncbi_taxonomy and ncbi_taxonomy != 'nan':
                        species_from_tax = extract_species_from_ncbi_taxonomy(ncbi_taxonomy)
                        if species_from_tax:
                            # Check if this species matches our TSV list
                            species_norm = normalize_species_name(species_from_tax).lower()
                            normalized_tsv_species = {normalize_species_name(s).lower(): s for s in tsv_species_list}
                            
                            if species_norm in normalized_tsv_species:
                                return normalized_tsv_species[species_norm]
                            
                            # Try partial matching
                            for norm_tsv_name, orig_tsv_name in normalized_tsv_species.items():
                                if species_norm in norm_tsv_name or norm_tsv_name in species_norm:
                                    return orig_tsv_name
                            
                            # Try Candidatus special case
                            candidatus_match = handle_candidatus_special_case(species_from_tax, tsv_species_list)
                            if candidatus_match:
                                return candidatus_match
                    
                    # If NCBI taxonomy didn't work, try organism name columns
                    for org_col in mapping_df.columns:
                        org_col_lower = org_col.lower()
                        if any(keyword in org_col_lower for keyword in ['organism', 'species', 'name']):
                            if pd.notna(row[org_col]) and str(row[org_col]).strip():
                                species_candidate = str(row[org_col]).strip()
                                if len(species_candidate.split()) >= 2:  # Proper binomial name
                                    species_norm = normalize_species_name(species_candidate).lower()
                                    normalized_tsv_species = {normalize_species_name(s).lower(): s for s in tsv_species_list}
                                    
                                    if species_norm in normalized_tsv_species:
                                        return normalized_tsv_species[species_norm]
                                    
                                    # Try partial matching
                                    for norm_tsv_name, orig_tsv_name in normalized_tsv_species.items():
                                        if species_norm in norm_tsv_name or norm_tsv_name in species_norm:
                                            return orig_tsv_name
                                    
                                    # Try Candidatus special case
                                    candidatus_match = handle_candidatus_special_case(species_candidate, tsv_species_list)
                                    if candidatus_match:
                                        return candidatus_match
    
    return None

def get_latest_version_files(fasta_dir):
    """Get the latest version of each protein file (prioritize r01_a01 > r01 > r00)."""
    # Find all nogaps files including a01 versions
    pattern1 = glob.glob(os.path.join(fasta_dir, "*_nogaps_r01_a01.fa"))  # Add this pattern
    pattern2 = glob.glob(os.path.join(fasta_dir, "*_nogaps_r01.fa"))
    pattern3 = glob.glob(os.path.join(fasta_dir, "*_nogaps_r00.fa"))
    
    all_files = pattern1 + pattern2 + pattern3
    
    # Group by protein name
    protein_files = defaultdict(list)
    for file_path in all_files:
        filename = os.path.basename(file_path)
        
        # Extract protein name more carefully
        if '_nogaps_r01_a01' in filename:
            protein_name = filename.replace('_nogaps_r01_a01.fa', '').split('_')[0]
        elif '_nogaps_r01' in filename:
            protein_name = filename.replace('_nogaps_r01.fa', '').split('_')[0]
        elif '_nogaps_r00' in filename:
            protein_name = filename.replace('_nogaps_r00.fa', '').split('_')[0]
        else:
            protein_name = filename.split('_')[0]
        
        protein_files[protein_name].append(file_path)
    
    # Select the latest version for each protein
    latest_files = []
    for protein, files in protein_files.items():
        # Sort to prioritize r01_a01 > r01 > r00
        def sort_key(file_path):
            filename = os.path.basename(file_path)
            if '_nogaps_r01_a01' in filename:
                return 3  # Highest priority
            elif '_nogaps_r01' in filename:
                return 2  # Medium priority
            elif '_nogaps_r00' in filename:
                return 1  # Lowest priority
            else:
                return 0  # Unknown format
        
        files.sort(key=sort_key, reverse=True)
        latest_files.append(files[0])
        print(f"Protein {protein}: Using {os.path.basename(files[0])}")
        if len(files) > 1:
            print(f"  (Skipping {[os.path.basename(f) for f in files[1:]]})")
    
    return latest_files

def get_combined_protein_files(combined_dir):
    """Get combined protein files from the combined_tree directory."""
    print(f"Looking for combined files in: {combined_dir}")
    
    # Find combined files with r01 and r00 versions
    pattern1 = glob.glob(os.path.join(combined_dir, "*_r01_a01.fa"))
    pattern2 = glob.glob(os.path.join(combined_dir, "*_r01.fa"))
    pattern3 = glob.glob(os.path.join(combined_dir, "*_r00.fa"))
    
    print(f"Pattern 1 (*_r01_a01.fa): {[os.path.basename(f) for f in pattern1]}")
    print(f"Pattern 2 (*_r01.fa): {[os.path.basename(f) for f in pattern2]}")
    print(f"Pattern 3 (*_r00.fa): {[os.path.basename(f) for f in pattern3]}")
    
    all_files = pattern1 + pattern2 + pattern3
    
    print(f"All combined files found: {[os.path.basename(f) for f in all_files]}")
    
    if not all_files:
        print("No combined files found!")
        return []
    
    # Group by combined protein name
    combined_protein_files = defaultdict(list)
    for file_path in all_files:
        filename = os.path.basename(file_path)
        
        # Extract combined name more carefully
        if '_r01_a01' in filename:
            combined_name = filename.replace('_r01_a01.fa', '')
        elif '_r01' in filename:
            combined_name = filename.replace('_r01.fa', '')
        elif '_r00' in filename:
            combined_name = filename.replace('_r00.fa', '')
        else:
            combined_name = filename.split('_')[0]  # fallback
        
        combined_protein_files[combined_name].append(file_path)
        print(f"  Grouped {filename} under protein name: '{combined_name}'")
    
    print(f"Protein groups created: {list(combined_protein_files.keys())}")
    
    # Select the latest version for each combined protein
    latest_combined_files = []
    for combined_protein, files in combined_protein_files.items():
        print(f"\nProcessing combined protein: '{combined_protein}'")
        print(f"Available files: {[os.path.basename(f) for f in files]}")
        
        # Sort to prioritize: r01_a01 > r01 > r00
        def sort_key(file_path):
            filename = os.path.basename(file_path)
            if '_r01_a01' in filename:
                return 3  # Highest priority
            elif '_r01' in filename and '_a01' not in filename:
                return 2  # Medium priority
            elif '_r00' in filename:
                return 1  # Lowest priority
            else:
                return 0  # Unknown format
        
        files_with_priority = [(os.path.basename(f), sort_key(f)) for f in files]
        print(f"Files with priorities: {files_with_priority}")
        
        files.sort(key=sort_key, reverse=True)
        selected_file = files[0]
        latest_combined_files.append(selected_file)
        print(f"SELECTED: {os.path.basename(selected_file)}")
        
        if len(files) > 1:
            print(f"SKIPPED: {[os.path.basename(f) for f in files[1:]]}")
    
    print(f"\nFINAL COMBINED FILES TO PROCESS: {[os.path.basename(f) for f in latest_combined_files]}")
    return latest_combined_files


    """Get combined protein files from the combined_tree directory."""
    # Find combined files with r01 and r00 versions
    pattern1 = glob.glob(os.path.join(combined_dir, "*_r01_a01.fa"))
    pattern2 = glob.glob(os.path.join(combined_dir, "*_r01.fa"))
    pattern3 = glob.glob(os.path.join(combined_dir, "*_r00.fa"))
    
    all_files = pattern1 + pattern2 + pattern3
    
    print(f"Found combined files:")
    for f in all_files:
        print(f"  {os.path.basename(f)}")
    
    # Group by combined protein name
    combined_protein_files = defaultdict(list)
    for file_path in all_files:
        filename = os.path.basename(file_path)
        
        # Extract combined name more carefully
        # Handle cases like: Tuc1_r01_a01.fa, MnmA-MTU1_r01.fa
        if '_r01_a01' in filename:
            combined_name = filename.replace('_r01_a01.fa', '')
        elif '_r01' in filename:
            combined_name = filename.replace('_r01.fa', '')
        elif '_r00' in filename:
            combined_name = filename.replace('_r00.fa', '')
        else:
            combined_name = filename.split('_')[0]  # fallback
        
        combined_protein_files[combined_name].append(file_path)
        print(f"  Grouped {filename} under protein name: {combined_name}")
    
    # Select the latest version for each combined protein
    latest_combined_files = []
    for combined_protein, files in combined_protein_files.items():
        print(f"\nProcessing combined protein: {combined_protein}")
        print(f"Available files: {[os.path.basename(f) for f in files]}")
        
        # Sort to prioritize: r01_a01 > r01 > r00
        def sort_key(file_path):
            filename = os.path.basename(file_path)
            if '_r01_a01' in filename:
                return 3  # Highest priority
            elif '_r01' in filename and '_a01' not in filename:
                return 2  # Medium priority
            elif '_r00' in filename:
                return 1  # Lowest priority
            else:
                return 0  # Unknown format
        
        files_with_priority = [(f, sort_key(f)) for f in files]
        print(f"Files with priorities: {[(os.path.basename(f), p) for f, p in files_with_priority]}")
        
        files.sort(key=sort_key, reverse=True)
        selected_file = files[0]
        latest_combined_files.append(selected_file)
        print(f"Selected: {os.path.basename(selected_file)}")
        
        if len(files) > 1:
            print(f"Skipped: {[os.path.basename(f) for f in files[1:]]}")
    
    return latest_combined_files

    """Get combined protein files from the combined_tree directory."""
    # Find combined files with r01 and r00 versions
    pattern1 = glob.glob(os.path.join(combined_dir, "*_r01_a01.fa"))
    pattern2 = glob.glob(os.path.join(combined_dir, "*_r01.fa"))
    pattern3 = glob.glob(os.path.join(combined_dir, "*_r00.fa"))
    
    all_files = pattern1 + pattern2 + pattern3
    
    # Group by combined protein name
    combined_protein_files = defaultdict(list)
    for file_path in all_files:
        filename = os.path.basename(file_path)
        combined_name = filename.split('_')[0]  # e.g., "MnmA-MTU1"
        combined_protein_files[combined_name].append(file_path)
    
    # Select the latest version for each combined protein
    latest_combined_files = []
    for combined_protein, files in combined_protein_files.items():
        # Sort to prioritize: r01_a01 > r01 > r00
        def sort_key(file_path):
            filename = os.path.basename(file_path)
            if '_r01_a01' in filename:
                return 3  # Highest priority
            elif '_r01' in filename:
                return 2  # Medium priority
            elif '_r00' in filename:
                return 1  # Lowest priority
            else:
                return 0  # Unknown format
        
        files.sort(key=sort_key, reverse=True)
        latest_combined_files.append(files[0])
        print(f"Combined protein {combined_protein}: Using {os.path.basename(files[0])}")
        if len(files) > 1:
            print(f"  (Skipping {[os.path.basename(f) for f in files[1:]]})")
    
    return latest_combined_files

def get_combined_protein_components(combined_name):
    """Extract individual protein components from combined protein name."""
    # Split by dash to get individual proteins
    components = combined_name.split('-')
    return components

def analyze_fasta_files(fasta_dir, combined_dir, tsv_species_list, protein_genome_mapping, protein_id_mapping, mapping_df):
    """Analyze all fasta files in the directory and return protein presence data."""
    # Get single protein files
    single_fasta_files = get_latest_version_files(fasta_dir)
    
    # Get combined protein files
    combined_fasta_files = get_combined_protein_files(combined_dir)
    
    print(f"Looking in single protein directory: {fasta_dir}")
    print(f"Found {len(single_fasta_files)} single protein files to process")
    print(f"Looking in combined protein directory: {combined_dir}")
    print(f"Found {len(combined_fasta_files)} combined protein files to process")
    print(f"Combined files: {[os.path.basename(f) for f in combined_fasta_files]}")
    
    # ... rest of the function stays the same
    # Create normalized species lookup
    normalized_tsv_species = {}
    for species in tsv_species_list:
        norm_name = normalize_species_name(species).lower()
        normalized_tsv_species[norm_name] = species
    
    print(f"Sample TSV species: {list(tsv_species_list)[:5]}")
    print(f"Sample normalized species: {list(normalized_tsv_species.keys())[:5]}")
    
    # Dictionary to store: {species_name: {protein_name: count}}
    species_protein_data = defaultdict(lambda: defaultdict(int))
    
    # Get all protein names from filenames
    single_protein_names = set()
    combined_protein_names = set()
    combined_protein_components = {}  # Maps combined protein -> list of components
    
    # Statistics
    total_matched_by_direct = 0
    total_matched_by_mapping = 0
    total_matched_by_protein_id = 0
    total_matched_by_direct_lookup = 0
    total_matched_by_candidatus = 0
    total_unmatched = 0
    
    # Process combined proteins first
    print("\n=== PROCESSING COMBINED PROTEINS ===")
    for fasta_file in combined_fasta_files:
        filename = os.path.basename(fasta_file)
        combined_protein_name = filename.split('_')[0]
        combined_protein_names.add(combined_protein_name)
        
        # Extract components
        components = get_combined_protein_components(combined_protein_name)
        combined_protein_components[combined_protein_name] = components
        print(f"Combined protein {combined_protein_name} contains: {components}")
        
        # Process this combined protein file
        species_counts = process_single_fasta_file(fasta_file, combined_protein_name, 
                                                   normalized_tsv_species, tsv_species_list,
                                                   protein_genome_mapping, protein_id_mapping, mapping_df)
        
        # Store the counts
        for species, count in species_counts.items():
            species_protein_data[species][combined_protein_name] = count
    
    # Process single proteins, but exclude those that are part of combined proteins
    print("\n=== PROCESSING SINGLE PROTEINS ===")
    all_combined_components = set()
    for components in combined_protein_components.values():
        all_combined_components.update(components)
    
    for fasta_file in single_fasta_files:
        filename = os.path.basename(fasta_file)
        protein_name = filename.split('_')[0]
        
        # Skip if this protein is part of a combined protein
        if protein_name in all_combined_components:
            print(f"Skipping {protein_name} - part of combined protein")
            continue
        
        single_protein_names.add(protein_name)
        
        # Process this single protein file
        species_counts = process_single_fasta_file(fasta_file, protein_name, 
                                                   normalized_tsv_species, tsv_species_list,
                                                   protein_genome_mapping, protein_id_mapping, mapping_df)
        
        # Store the counts
        for species, count in species_counts.items():
            species_protein_data[species][protein_name] = count
    
    # Combine all protein names
    all_protein_names = sorted(list(single_protein_names) + list(combined_protein_names))
    
    print(f"\nFound {len(single_protein_names)} single proteins: {sorted(single_protein_names)}")
    print(f"Found {len(combined_protein_names)} combined proteins: {sorted(combined_protein_names)}")
    print(f"Excluded single proteins (part of combined): {sorted(all_combined_components)}")
    
    return species_protein_data, all_protein_names

def process_single_fasta_file(fasta_file, protein_name, normalized_tsv_species, tsv_species_list,
                              protein_genome_mapping, protein_id_mapping, mapping_df):
    """Process a single FASTA file and return species counts."""
    print(f"\nProcessing {os.path.basename(fasta_file)} -> protein: {protein_name}")
    
    # Count sequences per species in this fasta file
    species_counts = defaultdict(int)
    species_sequences = defaultdict(list)  # Track actual sequences for debugging
    found_species = set()
    unmatched_examples = []
    matched_by_mapping = 0
    matched_by_direct = 0
    matched_by_protein_id = 0
    matched_by_direct_lookup = 0
    matched_by_candidatus = 0
    
    try:
        with open(fasta_file, 'r') as f:
            for line_num, line in enumerate(f):
                if line.startswith('>'):
                    # Extract species name from header
                    species_raw = extract_species_from_header(line)
                    species_norm = normalize_species_name(species_raw).lower()
                    
                    # Also extract protein ID
                    protein_id = extract_protein_id_from_header(line)
                    
                    # Try to match with TSV species
                    matched_species = None
                    
                    # Strategy 1: Direct matching with TSV species
                    if species_norm in normalized_tsv_species:
                        matched_species = normalized_tsv_species[species_norm]
                        matched_by_direct += 1
                    else:
                        # Try partial matching
                        for norm_tsv_name, orig_tsv_name in normalized_tsv_species.items():
                            if species_norm in norm_tsv_name or norm_tsv_name in species_norm:
                                matched_species = orig_tsv_name
                                matched_by_direct += 1
                                break
                    
                    # Strategy 2: Try Candidatus special case matching
                    if not matched_species and species_raw:
                        candidatus_match = handle_candidatus_special_case(species_raw, tsv_species_list)
                        if candidatus_match:
                            matched_species = candidatus_match
                            matched_by_candidatus += 1
                            print(f"    Matched by Candidatus special case: {species_raw} -> {matched_species}")
                    
                    # Strategy 3: If not found, try protein ID mapping (pre-built dictionary)
                    if not matched_species and protein_id and protein_id_mapping:
                        if protein_id in protein_id_mapping:
                            mapped_species = protein_id_mapping[protein_id]
                            mapped_norm = normalize_species_name(mapped_species).lower()
                            if mapped_norm in normalized_tsv_species:
                                matched_species = normalized_tsv_species[mapped_norm]
                                matched_by_protein_id += 1
                                print(f"    Matched by protein ID dict {protein_id}: {mapped_species}")
                            else:
                                # Try Candidatus special case on mapped species
                                candidatus_match = handle_candidatus_special_case(mapped_species, tsv_species_list)
                                if candidatus_match:
                                    matched_species = candidatus_match
                                    matched_by_candidatus += 1
                                    print(f"    Matched by protein ID + Candidatus {protein_id}: {mapped_species} -> {matched_species}")
                    
                    # Strategy 4: If still not found, try direct lookup in mapping dataframe
                    if not matched_species and protein_id:
                        direct_lookup_result = lookup_protein_id_in_mapping(protein_id, mapping_df, tsv_species_list)
                        if direct_lookup_result:
                            matched_species = direct_lookup_result
                            matched_by_direct_lookup += 1
                            print(f"    Matched by direct lookup {protein_id}: {matched_species}")
                    
                    # Strategy 5: If still not found, try protein genome mapping
                    if not matched_species and protein_genome_mapping:
                        # Try exact match first
                        if species_norm in protein_genome_mapping:
                            mapped_species = protein_genome_mapping[species_norm]
                            mapped_norm = normalize_species_name(mapped_species).lower()
                            if mapped_norm in normalized_tsv_species:
                                matched_species = normalized_tsv_species[mapped_norm]
                                matched_by_mapping += 1
                            else:
                                # Try Candidatus special case on mapped species
                                candidatus_match = handle_candidatus_special_case(mapped_species, tsv_species_list)
                                if candidatus_match:
                                    matched_species = candidatus_match
                                    matched_by_candidatus += 1
                                    print(f"    Matched by genome mapping + Candidatus: {mapped_species} -> {matched_species}")
                        else:
                            # Try partial matching with mapping keys
                            for mapping_key, mapping_species in protein_genome_mapping.items():
                                if species_norm in mapping_key or mapping_key in species_norm:
                                    mapped_norm = normalize_species_name(mapping_species).lower()
                                    if mapped_norm in normalized_tsv_species:
                                        matched_species = normalized_tsv_species[mapped_norm]
                                        matched_by_mapping += 1
                                        break
                                    else:
                                        # Try Candidatus special case on mapped species
                                        candidatus_match = handle_candidatus_special_case(mapping_species, tsv_species_list)
                                        if candidatus_match:
                                            matched_species = candidatus_match
                                            matched_by_candidatus += 1
                                            print(f"    Matched by genome mapping + Candidatus: {mapping_species} -> {matched_species}")
                                            break
                    
                    if matched_species:
                        species_counts[matched_species] += 1
                        species_sequences[matched_species].append(line.strip()[:80] + "...")  # Store header for debugging
                        found_species.add(species_raw)
                    else:
                        # Store unmatched examples (only first 10)
                        if len(unmatched_examples) < 10 and (species_raw or protein_id):
                            unmatched_examples.append({
                                'header': line.strip()[:100] + "...",
                                'extracted': species_raw,
                                'normalized': species_norm,
                                'protein_id': protein_id
                            })
                        
    except Exception as e:
        print(f"Error reading file {fasta_file}: {e}")
        return defaultdict(int)
    
    print(f"  Found {len(species_counts)} matched species with {sum(species_counts.values())} total sequences")
    print(f"  Matched by direct: {matched_by_direct}, by Candidatus: {matched_by_candidatus}")
    print(f"  Matched by protein ID dict: {matched_by_protein_id}, by direct lookup: {matched_by_direct_lookup}")
    print(f"  Matched by genome mapping: {matched_by_mapping}")
    print(f"  Unique species found: {list(found_species)[:5]}...")
    
    # Show species with high counts for debugging
    high_count_species = [(species, count) for species, count in species_counts.items() if count > 3]
    if high_count_species:
        print(f"  HIGH COUNT SPECIES (>3 sequences):")
        for species, count in sorted(high_count_species, key=lambda x: x[1], reverse=True)[:3]:
            print(f"    {species}: {count} sequences")
            for seq_header in species_sequences[species][:3]:  # Show first 3 headers
                print(f"      {seq_header}")
            if len(species_sequences[species]) > 3:
                print(f"      ... and {len(species_sequences[species]) - 3} more")
    
    # Show unmatched examples
    if unmatched_examples:
        print(f"  STILL UNMATCHED examples ({len(unmatched_examples)}):")
        for ex in unmatched_examples:
            print(f"    Header: {ex['header']}")
            print(f"    Extracted: '{ex['extracted']}' -> Normalized: '{ex['normalized']}'")
            if ex['protein_id']:
                print(f"    Protein ID: {ex['protein_id']}")
            print()
    
    # Convert set to dict for compatibility (each species present = 1)
    return species_counts

def debug_species_matching(species_name, tsv_species_list):
    """Debug function to see how a specific species is being matched."""
    print(f"\n=== DEBUGGING SPECIES MATCHING FOR: {species_name} ===")
    
    # Check normalization
    normalized = normalize_species_name(species_name).lower()
    print(f"Normalized: '{normalized}'")
    
    # Create normalized TSV lookup
    normalized_tsv_species = {}
    for species in tsv_species_list:
        norm_name = normalize_species_name(species).lower()
        normalized_tsv_species[norm_name] = species
    
    # Check if it matches directly
    if normalized in normalized_tsv_species:
        print(f"DIRECT MATCH: {normalized_tsv_species[normalized]}")
        return normalized_tsv_species[normalized]
    
    # Check partial matches
    partial_matches = []
    for norm_tsv_name, orig_tsv_name in normalized_tsv_species.items():
        if normalized in norm_tsv_name or norm_tsv_name in normalized:
            partial_matches.append((norm_tsv_name, orig_tsv_name))
    
    if partial_matches:
        print(f"PARTIAL MATCHES: {partial_matches}")
        return partial_matches[0][1]
    
    print("NO MATCH FOUND")
    return None

def main():
    # Read the original TSV file
    input_file = sys.argv[1]
    fasta_dir = sys.argv[2]
    combined_dir = sys.argv[3]
    mapping_file = sys.argv[4]
    
    # Load protein genome mapping
    print("Loading protein genome mapping...")
    protein_genome_mapping, protein_id_mapping, mapping_df = load_protein_genome_mapping(mapping_file)
    
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', header=None, 
                     names=['Domain', 'Phylum', 'Order', 'Family', 'Species'])
    
    print(f"Loaded {len(df)} species from TSV file")
    
    # Analyze fasta files (both single and combined)
    print("Analyzing fasta files...")
    species_protein_data, protein_names = analyze_fasta_files(fasta_dir, combined_dir, df['Species'].tolist(), 
                                                              protein_genome_mapping, protein_id_mapping, mapping_df)
    
    print(f"\nFound {len(protein_names)} total proteins: {protein_names}")
    print(f"Found data for {len(species_protein_data)} species")
    
    if len(protein_names) == 0:
        print("No proteins found. Exiting.")
        return
    
    # Add protein columns to dataframe
    for protein in protein_names:
        protein_column = []
        
        for _, row in df.iterrows():
            species = row['Species']
            
            # Check if species has this protein
            if species in species_protein_data and protein in species_protein_data[species]:
                count = species_protein_data[species][protein]
                protein_column.append(count)
            else:
                protein_column.append('FALSE')
        
        df[protein] = protein_column
    
    # Save the updated dataframe
    output_file = sys.argv[5]
    df.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"\nAnalysis complete! Output saved to: {output_file}")
    print(f"Added {len(protein_names)} protein columns to {len(df)} species")
    
    # Create and save statistics
    statistics_data = []
    for protein in protein_names:
        present_count = sum(1 for val in df[protein] if val != 'FALSE')
        total_copies = sum(val for val in df[protein] if val != 'FALSE')
        presence_percentage = (present_count / len(df)) * 100
        avg_copies_per_species = total_copies / present_count if present_count > 0 else 0
        
        statistics_data.append({
            'Protein': protein,
            'Present_in_species': present_count,
            'Total_species': len(df),
            'Presence_percentage': round(presence_percentage, 2),
            'Total_copies': total_copies,
            'Avg_copies_per_species': round(avg_copies_per_species, 2)
        })
        
        # Also print to console
        print(f"{protein}: Present in {present_count}/{len(df)} species ({presence_percentage:.2f}%), Total copies: {total_copies}, Avg: {avg_copies_per_species:.2f}")
    
    # Save statistics to file
    stats_df = pd.DataFrame(statistics_data)
    stats_file = sys.argv[6] 
    stats_df.to_csv(stats_file, sep='\t', index=False, header=True)
    
    print(f"\nProtein presence statistics saved to: {stats_file}")

if __name__ == "__main__":
    main()