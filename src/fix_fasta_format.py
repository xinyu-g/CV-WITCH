#!/usr/bin/env python3

import os
import sys
import re

def fix_fasta_format(input_file, output_file):
    """
    Fix FASTA file formatting issues:
    1. Ensure proper header/sequence separation
    2. Remove divider entries
    3. Format sequences with proper line breaks
    """
    with open(input_file, 'r') as f:
        content = f.read().strip()
    
    # Split into entries (assuming entries are separated by '>')
    entries = content.split('>')
    # Remove empty entries
    entries = [e.strip() for e in entries if e.strip()]
    
    with open(output_file, 'w') as f:
        for entry in entries:
            # Skip divider entries
            if '____DIVIDER' in entry:
                continue
                
            # Split into header and sequence
            lines = entry.split('\n')
            header = lines[0]
            # Join all remaining lines and remove any whitespace/newlines
            sequence = ''.join(lines[1:]).replace(' ', '').replace('\n', '')
            
            # Write header
            f.write(f">{header}\n")
            
            # Write sequence in chunks of 60 characters
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python fix_fasta_format.py <input_fasta> <output_fasta>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist")
        sys.exit(1)
        
    try:
        fix_fasta_format(input_file, output_file)
        print(f"Successfully reformatted {input_file} to {output_file}")
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 