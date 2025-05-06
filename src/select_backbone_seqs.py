#!/usr/bin/env python3

import sys
from Bio import SeqIO
import numpy as np
from pathlib import Path

def get_sequence_lengths(fasta_file):
    """Get sequence lengths and their IDs from FASTA file"""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((record.id, len(record.seq)))
    return sequences

def select_backbone_sequences(fasta_file, output_file, percentage=0.1):
    """
    Select sequences closest to median length as backbone sequences
    
    Args:
        fasta_file: Input FASTA file
        output_file: Output FASTA file for backbone sequences
        percentage: Percentage of sequences to select (default: 0.1 for 10%)
    """
    # Get sequence lengths
    sequences = get_sequence_lengths(fasta_file)
    seq_ids, lengths = zip(*sequences)
    
    # Calculate median length
    median_length = np.median(lengths)
    
    # Calculate distance from median for each sequence
    distances = [(seq_id, abs(length - median_length)) 
                for seq_id, length in sequences]
    
    # Sort by distance from median
    distances.sort(key=lambda x: x[1])
    
    # Select top percentage closest to median
    n_select = max(1, int(len(sequences) * percentage))
    selected_ids = set(x[0] for x in distances[:n_select])
    
    # Write selected sequences to output file
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in selected_ids:
                SeqIO.write(record, out_f, "fasta")
    
    print(f"Selected {n_select} sequences for backbone")
    return n_select

def main():
    if len(sys.argv) != 3:
        print("Usage: python select_backbone_seqs.py <input_fasta> <output_fasta>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    n_selected = select_backbone_sequences(input_file, output_file)
    print(f"Wrote {n_selected} backbone sequences to {output_file}")

if __name__ == "__main__":
    main() 