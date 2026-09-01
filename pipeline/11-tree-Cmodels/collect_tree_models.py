#!/usr/bin/env python3
"""
Collect tree models from IQ-TREE log files.
Parses _Cmod.log files from both regular and combined tree directories.
"""

import os
import re
import glob
import pandas as pd

def parse_iqtree_log(log_file):
    """
    Parse IQ-TREE log file to extract model information.
    Returns dict with model details or None if parsing fails.
    """
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Extract protein name from filename
        basename = os.path.basename(log_file)
        if basename.endswith('_Cmod.log'):
            protein_name = basename.replace('_Cmod.log', '')
        else:
            protein_name = basename.replace('.log', '')
        
        # Look for model information
        best_model = None
        log_likelihood = None
        aic_score = None
        bic_score = None
        
        # Try different patterns for model information
        
        # Pattern 1: Best-fit model (from model selection)
        best_model_pattern = r'Best-fit model:\s+(\S+)'
        model_match = re.search(best_model_pattern, content)
        if model_match:
            best_model = model_match.group(1)
        
        # Pattern 2: Model from command line (fixed model)
        if not best_model:
            command_pattern = r'iqtree.*?-m\s+(\S+)'
            command_match = re.search(command_pattern, content)
            if command_match:
                best_model = command_match.group(1)
        
        # Pattern 3: Model name from analysis section
        if not best_model:
            analysis_pattern = r'Model of evolution:\s+(\S+)'
            analysis_match = re.search(analysis_pattern, content)
            if analysis_match:
                best_model = analysis_match.group(1)
        
        # Extract log likelihood - try multiple patterns
        ll_patterns = [
            r'Log-likelihood of the tree:\s+([-\d\.]+)',
            r'Loglikelihood:\s+([-\d\.]+)',
            r'Log likelihood:\s+([-\d\.]+)'
        ]
        
        for pattern in ll_patterns:
            ll_match = re.search(pattern, content)
            if ll_match:
                log_likelihood = float(ll_match.group(1))
                break
        
        # Extract AIC
        aic_pattern = r'Akaike information criterion \(AIC\) score:\s+([\d\.]+)'
        aic_match = re.search(aic_pattern, content)
        if aic_match:
            aic_score = float(aic_match.group(1))
        
        # Extract BIC
        bic_pattern = r'Bayesian information criterion \(BIC\) score:\s+([\d\.]+)'
        bic_match = re.search(bic_pattern, content)
        if bic_match:
            bic_score = float(bic_match.group(1))
        
        # Check if this is a combined tree
        is_combined = 'combined_' in protein_name or '-' in protein_name
        
        return {
            'protein': protein_name,
            'type': 'combined' if is_combined else 'single',
            'best_model': best_model,
            'log_likelihood': log_likelihood,
            'aic_score': aic_score,
            'bic_score': bic_score,
            'log_file': log_file
        }
        
    except Exception as e:
        print(f"Error parsing {log_file}: {e}")
        return None

def collect_tree_models():
    """
    Collect tree models from all _Cmod.log files in specified directories.
    """
    # Define search paths
    search_paths = [
        '/home/begum/PycharmProjects/PythonProject/Jan_2026_tRNA/output/iqtree/combined_tree',
        '/home/begum/PycharmProjects/PythonProject/Jan_2026_tRNA/output/iqtree'
    ]
    
    all_results = []
    
    for search_path in search_paths:
        if not os.path.exists(search_path):
            print(f"Warning: Path does not exist: {search_path}")
            continue
        
        # Find all _Cmod.log files
        log_pattern = os.path.join(search_path, '*_Cmod.log')
        log_files = glob.glob(log_pattern)
        
        print(f"Found {len(log_files)} _Cmod.log files in {search_path}")
        
        for log_file in log_files:
            result = parse_iqtree_log(log_file)
            if result:
                all_results.append(result)
                print(f"  Parsed: {result['protein']} -> {result['best_model']}")
    
    return all_results

def save_results(results, output_file):
    """
    Save results to TSV file and print summary.
    """
    if not results:
        print("No results to save!")
        return
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Sort by type (combined first) then by protein name
    df = df.sort_values(['type', 'protein'])
    
    # Save to TSV
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Print summary
    print(f"\nSummary:")
    print(f"  Total proteins: {len(df)}")
    print(f"  Single proteins: {len(df[df['type'] == 'single'])}")
    print(f"  Combined proteins: {len(df[df['type'] == 'combined'])}")
    
    # Model frequency
    print(f"\nModel frequency:")
    model_counts = df['best_model'].value_counts()
    for model, count in model_counts.items():
        print(f"  {model}: {count}")
    
    # Show first few results
    print(f"\nFirst 10 results:")
    print(df[['protein', 'type', 'best_model', 'aic_score']].head(10).to_string(index=False))

def main():
    """Main function."""
    print("Collecting tree models from IQ-TREE log files...")
    
    # Collect results
    results = collect_tree_models()
    
    if results:
        # Output file
        output_file = 'tree_models_summary.tsv'
        save_results(results, output_file)
    else:
        print("No log files found or parsed successfully!")

if __name__ == "__main__":
    main()