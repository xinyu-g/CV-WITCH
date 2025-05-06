#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import argparse

def run_fasttree(alignment_file, output_tree, is_protein=False):
    """Run FastTree on the alignment file"""
    cmd = ['FastTree']
    if is_protein:
        cmd.append('-lg')  # Use LG model for protein sequences
    else:
        cmd.append('-nt')  # Use nucleotide model
    
    with open(alignment_file) as f_in, open(output_tree, 'w') as f_out:
        subprocess.run(cmd, stdin=f_in, stdout=f_out, check=True)

def calculate_sp_scores(test_alignment, ref_alignment):
    """
    Calculate Sum-of-Pairs scores (SPFN and SPFP) between test and reference alignments
    """
    def get_aligned_pairs(alignment):
        pairs = set()
        for i in range(len(alignment[0])):
            col = [seq[i] for seq in alignment]
            for j in range(len(col)):
                if col[j] != '-':
                    for k in range(j + 1, len(col)):
                        if col[k] != '-':
                            pairs.add((j, k, i))
        return pairs
    
    # Get aligned pairs from both alignments
    test_pairs = get_aligned_pairs(test_alignment)
    ref_pairs = get_aligned_pairs(ref_alignment)
    
    # Calculate metrics
    true_positives = len(test_pairs.intersection(ref_pairs))
    false_positives = len(test_pairs - ref_pairs)
    false_negatives = len(ref_pairs - test_pairs)
    
    # Calculate rates
    spfn = false_negatives / len(ref_pairs) if ref_pairs else 0
    spfp = false_positives / len(test_pairs) if test_pairs else 0
    
    return spfn, spfp

def calculate_fn_rate(test_tree_file, ref_tree_file):
    """
    Calculate False Negative (FN) rate by comparing trees
    Uses Robinson-Foulds distance normalized by tree size
    """
    # Read trees using FastTree's built-in comparison
    cmd = ['FastTree', '-compare', ref_tree_file, test_tree_file]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
    # Parse the output to get RF distance
    for line in result.stdout.split('\n'):
        if 'Robinson-Foulds' in line:
            rf_dist = float(line.split()[-1])
            break
    
    # Normalize by tree size (number of internal nodes)
    with open(test_tree_file) as f:
        tree_str = f.read()
        n_leaves = tree_str.count(',') + 1
        n_internal = n_leaves - 2  # Number of internal nodes in a binary tree
    
    fn_rate = rf_dist / (2 * n_internal)  # Normalize by maximum possible RF distance
    return fn_rate

def evaluate_alignment(test_aln_file, ref_aln_file, output_dir, is_protein=False):
    """
    Evaluate an alignment by comparing it to a reference alignment
    """
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read alignments
    test_aln = AlignIO.read(test_aln_file, "fasta")
    ref_aln = AlignIO.read(ref_aln_file, "fasta")
    
    # Calculate SP scores
    spfn, spfp = calculate_sp_scores(test_aln, ref_aln)
    
    # Generate trees
    test_tree = output_dir / "test_tree.nwk"
    ref_tree = output_dir / "ref_tree.nwk"
    
    run_fasttree(test_aln_file, test_tree, is_protein)
    run_fasttree(ref_aln_file, ref_tree, is_protein)
    
    # Calculate FN rate
    fn_rate = calculate_fn_rate(str(test_tree), str(ref_tree))
    
    return {
        'SPFN': spfn,
        'SPFP': spfp,
        'FN': fn_rate
    }

def main():
    parser = argparse.ArgumentParser(description='Evaluate multiple sequence alignments using FastTree')
    parser.add_argument('--test-dir', required=True, help='Directory containing test alignments')
    parser.add_argument('--ref-dir', required=True, help='Directory containing reference alignments')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    parser.add_argument('--protein', action='store_true', help='Use protein model (default: nucleotide)')
    args = parser.parse_args()
    
    test_dir = Path(args.test_dir)
    ref_dir = Path(args.ref_dir)
    output_dir = Path(args.output_dir)
    
    # Store results for each dataset
    results = {}
    
    # Process all alignment files
    for test_aln_file in test_dir.glob("**/*.fasta"):
        # Find corresponding reference alignment
        rel_path = test_aln_file.relative_to(test_dir)
        ref_aln_file = ref_dir / rel_path
        
        if not ref_aln_file.exists():
            print(f"Warning: No reference alignment found for {test_aln_file}")
            continue
        
        print(f"Processing {rel_path}...")
        
        # Create output directory for this alignment
        aln_output_dir = output_dir / rel_path.parent
        
        try:
            metrics = evaluate_alignment(
                test_aln_file,
                ref_aln_file,
                aln_output_dir,
                is_protein=args.protein
            )
            
            results[str(rel_path)] = metrics
            
            # Write individual results
            with open(aln_output_dir / f"{rel_path.stem}_metrics.txt", 'w') as f:
                for metric, value in metrics.items():
                    f.write(f"{metric}: {value:.4f}\n")
            
        except Exception as e:
            print(f"Error processing {rel_path}: {e}")
            continue
    
    # Write summary results
    with open(output_dir / "summary_metrics.txt", 'w') as f:
        f.write("Alignment\tSPFN\tSPFP\tFN\n")
        for aln_path, metrics in results.items():
            f.write(f"{aln_path}\t{metrics['SPFN']:.4f}\t{metrics['SPFP']:.4f}\t{metrics['FN']:.4f}\n")
    
    # Calculate and write average metrics
    avg_metrics = {
        metric: np.mean([r[metric] for r in results.values()])
        for metric in ['SPFN', 'SPFP', 'FN']
    }
    
    with open(output_dir / "average_metrics.txt", 'w') as f:
        for metric, value in avg_metrics.items():
            f.write(f"Average {metric}: {value:.4f}\n")

if __name__ == "__main__":
    main() 