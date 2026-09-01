#!/usr/bin/env python3
import re
import sys

def add_logs_to_snakefile(filename):
    with open(filename, 'r') as f:
        content = f.read()
    
    # Pattern to match rules without log directive
    # Matches: rule name: ... output: ... shell/run:
    rule_pattern = r'(rule\s+(\w+):\s*\n(?:.*\n)*?)(    shell:)'
    
    def add_log(match):
        rule_block = match.group(1)
        rule_name = match.group(2)
        shell_line = match.group(3)
        
        # Skip if already has log
        if 'log:' in rule_block:
            return match.group(0)
        
        # Extract wildcard from output line
        output_match = re.search(r'output:.*?{(\w+)}', rule_block)
        wildcard = output_match.group(1) if output_match else 'output'
        
        # Add log directive before shell
        log_line = f'    log:\n        logs + "{rule_name}_{{{wildcard}}}.log"\n'
        
        # Update shell command to redirect to log
        updated_shell = rule_block + log_line + shell_line
        
        return updated_shell
    
    # Replace all rules
    content = re.sub(rule_pattern, add_log, content)
    
    # Now update shell commands to append | tee {log}
    # For simple shell commands
    content = re.sub(
        r'(    log:\n        logs.*\n    shell:\n        ")(.*?)(")',
        lambda m: m.group(1) + m.group(2) + ' 2>&1 | tee {log}' + m.group(3),
        content
    )
    
    # For triple-quoted shell commands
    content = re.sub(
        r'(    log:\n        logs.*\n    shell:\n        """)(.*?)(\n        """)',
        lambda m: m.group(1) + m.group(2) + ' 2>&1 | tee {log}' + m.group(3),
        content,
        flags=re.DOTALL
    )
    
    with open(filename, 'w') as f:
        f.write(content)
    
    print(f"Updated {filename}")

if __name__ == '__main__':
    for filename in ['blast-all-in-bact-db.smk', 'blast-all-in-euk-db.smk']:
        add_logs_to_snakefile(filename)
