#!/usr/bin/env python3
"""
Update taxonomic information for Archaea and Bacteria based on merged_meta_data_bact_arc.tsv
The order is determined by dataset_taxa_order_manual.txt, but taxonomy comes from metadata.
"""

import pandas as pd
import sys
import os

def main():
    # Read the metadata table from command line argument or use relative path as fallback
    if len(sys.argv) > 1:
        metadata_file = sys.argv[1]
        input_file = sys.argv[2]
        output_file = sys.argv[3]
   
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    # Create a mapping from species to full taxonomy row
    # We'll store both the parsed taxonomy and the original metadata row
    species_to_taxonomy = {}
    
    for idx, row in metadata.iterrows():
        # Get NCBI taxonomy data and organism name
        ncbi_taxonomy = str(row.get('ncbi_taxonomy', ''))
        ncbi_organism = str(row.get('ncbi_organism_name', ''))
        
        # Skip if no NCBI taxonomy available
        if not ncbi_taxonomy or ncbi_taxonomy == 'nan':
            continue
            
        # Parse NCBI taxonomy string to extract components
        # NCBI format: d__;p__;c__;o__;f__;g__;s__ (similar to GTDB but from NCBI)
        parts = ncbi_taxonomy.split(';')
        
        # Extract domain, phylum, order, family, species
        domain = 'N/A'
        phylum = 'N/A'
        order = 'N/A'
        family = 'N/A'
        species_from_taxonomy = 'N/A'
        
        for part in parts:
            part = part.strip()
            if part.startswith('d__'):
                extracted = part.replace('d__', '').strip()
                domain = extracted if extracted else 'N/A'
            elif part.startswith('p__'):
                extracted = part.replace('p__', '').strip()
                phylum = extracted if extracted else 'N/A'
            elif part.startswith('o__'):
                extracted = part.replace('o__', '').strip()
                order = extracted if extracted else 'N/A'
            elif part.startswith('f__'):
                extracted = part.replace('f__', '').strip()
                family = extracted if extracted else 'N/A'
            elif part.startswith('s__'):
                extracted = part.replace('s__', '').strip()
                species_from_taxonomy = extracted if extracted else 'N/A'
        
        # Use NCBI organism name ONLY if species from taxonomy is empty/N/A
        if species_from_taxonomy == 'N/A' or not species_from_taxonomy:
            final_species_name = ncbi_organism if (ncbi_organism and ncbi_organism != 'nan') else 'N/A'
        else:
            final_species_name = species_from_taxonomy
        
        # Create the taxonomy entry
        tax_entry = {
            'domain': domain,
            'phylum': phylum,
            'order': order,
            'family': family,
            'species': final_species_name
        }
        
        # Store with multiple keys for better matching
        # 1. Full ncbi_organism_name (if available)
        if ncbi_organism and ncbi_organism != 'nan':
            species_to_taxonomy[ncbi_organism.lower()] = tax_entry
        
        # 2. Full species from taxonomy (if available)
        if final_species_name and final_species_name != 'N/A':
            species_to_taxonomy[final_species_name.lower()] = tax_entry
        
        # 3. Genus + species (first two words)
        if final_species_name and final_species_name != 'N/A':
            species_parts = final_species_name.split()
            if len(species_parts) >= 2:
                genus_species = f"{species_parts[0]} {species_parts[1]}"
                species_to_taxonomy[genus_species.lower()] = tax_entry
    
    # Read the species order file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    updated_lines = []
    updates_count = 0
    not_found = []
    
    for line in lines:
        line = line.rstrip('\n')
        parts = line.split('\t')
        
        if len(parts) < 5:
            updated_lines.append(line)
            continue
        
        domain = parts[0]
        species = parts[4]
        
        # For Eukarya, keep as is (replace underscores with spaces)
        if domain == 'Eukarya':
            species_clean = species.replace('_', ' ')
            updated_line = f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{species_clean}"
            updated_lines.append(updated_line)
            continue
        
        # For Archaea and Bacteria, lookup taxonomy from metadata
        if domain in ['Archaea', 'Bacteria']:
            matched = False
            
            # Try different matching strategies
            species_lower = species.lower()
            
            # 1. Exact match
            if species_lower in species_to_taxonomy:
                tax = species_to_taxonomy[species_lower]
                updated_line = f"{tax['domain']}\t{tax['phylum']}\t{tax['order']}\t{tax['family']}\t{tax['species']}"
                updated_lines.append(updated_line)
                updates_count += 1
                matched = True
            else:
                # 2. Try genus-only match (first word)
                species_words = species.split()
                if len(species_words) >= 1:
                    genus = species_words[0].lower()
                    
                    # Look for any entry starting with this genus
                    for key, tax in species_to_taxonomy.items():
                        if key.startswith(genus + ' '):
                            updated_line = f"{tax['domain']}\t{tax['phylum']}\t{tax['order']}\t{tax['family']}\t{tax['species']}"
                            updated_lines.append(updated_line)
                            updates_count += 1
                            matched = True
                            break
                
                # 3. Try partial match (in case of strain differences)
                if not matched and len(species_words) >= 2:
                    genus_species = f"{species_words[0]} {species_words[1]}".lower()
                    
                    for key, tax in species_to_taxonomy.items():
                        if genus_species in key or key.startswith(genus_species):
                            updated_line = f"{tax['domain']}\t{tax['phylum']}\t{tax['order']}\t{tax['family']}\t{tax['species']}"
                            updated_lines.append(updated_line)
                            updates_count += 1
                            matched = True
                            break
                
                # 4. Special case: similar genus names (e.g., Borreliella vs Borrelia)
                # Match by genus prefix (first 5-6 characters)
                if not matched and len(species_words) >= 1:
                    genus_prefix = genus[:6] if len(genus) > 6 else genus[:5]
                    
                    for key, tax in species_to_taxonomy.items():
                        key_genus = key.split()[0] if ' ' in key else key
                        if len(key_genus) > 4 and key_genus[:6].startswith(genus_prefix):
                            updated_line = f"{tax['domain']}\t{tax['phylum']}\t{tax['order']}\t{tax['family']}\t{tax['species']}"
                            updated_lines.append(updated_line)
                            updates_count += 1
                            matched = True
                            break
            
            if not matched:
                # Keep original line but mark fields that were empty as N/A
                tax_parts = parts.copy()
                for i in range(1, 4):  # phylum, order, family
                    if not tax_parts[i] or tax_parts[i].strip() == '' or tax_parts[i].strip() == 'False':
                        tax_parts[i] = 'N/A'
                updated_line = '\t'.join(tax_parts)
                updated_lines.append(updated_line)
                not_found.append(species)
        else:
            # Keep other domains unchanged
            updated_lines.append(line)
    
    # Write updated file
    with open(output_file, 'w') as f:
        for line in updated_lines:
            f.write(line + '\n')
    
    print(f"Updated {updates_count} entries")
    print(f"Total lines: {len(lines)}")
    if not_found:
        print(f"\nNot found in metadata ({len(not_found)} entries):")
        for sp in not_found[:10]:  # Show first 10
            print(f"  - {sp}")
        if len(not_found) > 10:
            print(f"  ... and {len(not_found) - 10} more")

if __name__ == '__main__':
    main()