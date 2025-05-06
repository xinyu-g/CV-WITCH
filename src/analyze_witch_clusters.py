#!/usr/bin/env python3

import os
import sys
from pathlib import Path
import json
from collections import defaultdict

def analyze_cluster_contents(cluster_dir):
    """Analyze contents of a single cluster directory"""
    info = {
        'sequences': [],
        'has_alignment': False,
        'has_hmm': False
    }
    
    # Check for sequence file
    seq_file = cluster_dir / "sequences.fasta"
    if seq_file.exists():
        with open(seq_file) as f:
            for line in f:
                if line.startswith('>'):
                    info['sequences'].append(line.strip()[1:])  # Remove '>'
    
    # Check for alignment and HMM
    info['has_alignment'] = (cluster_dir / "alignment.fasta").exists()
    info['has_hmm'] = (cluster_dir / "hmm.txt").exists()
    
    return info

def analyze_tree_structure(tree_file):
    """Analyze the decomposition tree structure"""
    if not tree_file.exists():
        return None
    
    tree_info = {
        'raw_tree': '',
        'relationships': defaultdict(list)
    }
    
    with open(tree_file) as f:
        tree_info['raw_tree'] = f.read()
        # Parse tree file to extract parent-child relationships
        # This depends on the exact format of decomposition.tree
    
    return tree_info

def analyze_weights(weights_file):
    """
    Analyze the weights file from WITCH output
    
    The weights file contains adjusted bitscore weights for each HMM-sequence pair,
    indicating how much each HMM's alignment should contribute to the final alignment
    for each sequence.
    """
    if not weights_file.exists():
        return None
        
    weights_info = {
        'hmm_weights': defaultdict(dict),  # HMM -> sequence -> weight
        'sequence_coverage': defaultdict(list),  # sequence -> list of HMMs
        'hmm_usage': defaultdict(int)  # HMM -> number of sequences it aligns
    }
    
    with open(weights_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                # Format: sequence_id hmm_id weight
                parts = line.strip().split()
                if len(parts) >= 3:
                    seq_id, hmm_id, weight = parts[0], parts[1], float(parts[2])
                    
                    # Store weight
                    weights_info['hmm_weights'][hmm_id][seq_id] = weight
                    
                    # Track which HMMs align each sequence
                    weights_info['sequence_coverage'][seq_id].append(hmm_id)
                    
                    # Count how many sequences each HMM aligns
                    weights_info['hmm_usage'][hmm_id] += 1
    
    return weights_info

def write_weights_analysis(weights_info, output_file):
    """Write detailed analysis of weights"""
    if not weights_info:
        return
        
    with open(output_file, 'w') as f:
        f.write("WITCH Weights Analysis\n")
        f.write("=====================\n\n")
        
        # Overall statistics
        f.write("Summary Statistics:\n")
        f.write("-----------------\n")
        f.write(f"Total HMMs: {len(weights_info['hmm_weights'])}\n")
        f.write(f"Total sequences: {len(weights_info['sequence_coverage'])}\n\n")
        
        # HMM usage analysis
        f.write("HMM Usage:\n")
        f.write("----------\n")
        for hmm_id, count in sorted(weights_info['hmm_usage'].items(), key=lambda x: x[1], reverse=True):
            f.write(f"HMM {hmm_id}: used for {count} sequences\n")
        f.write("\n")
        
        # Sequence coverage analysis
        f.write("Sequence Coverage:\n")
        f.write("-----------------\n")
        for seq_id, hmms in sorted(weights_info['sequence_coverage'].items()):
            f.write(f"Sequence {seq_id}:\n")
            f.write(f"  Aligned by {len(hmms)} HMMs\n")
            for hmm_id in hmms:
                weight = weights_info['hmm_weights'][hmm_id][seq_id]
                f.write(f"  - HMM {hmm_id}: weight = {weight:.4f}\n")
            f.write("\n")

def analyze_witch_output(witch_output_dir):
    """
    Analyze WITCH's cluster structure and merging process
    
    Args:
        witch_output_dir: Path to WITCH output directory
    """
    output_dir = Path(witch_output_dir)
    decomp_dir = output_dir / "tree_decomp"
    root_dir = decomp_dir / "root"
    weights_file = output_dir / "weights.txt"
    
    if not root_dir.exists():
        print(f"Error: Cannot find tree_decomp/root in {witch_output_dir}")
        return
    
    # Analyze clusters
    clusters = {}
    for cluster_dir in root_dir.glob("A_*"):
        clusters[cluster_dir.name] = analyze_cluster_contents(cluster_dir)
    
    # Analyze tree structure
    tree_info = analyze_tree_structure(decomp_dir / "decomposition.tree")
    
    # Analyze weights
    weights_info = analyze_weights(weights_file)
    
    # Generate reports
    report = {
        'summary': {
            'total_clusters': len(clusters),
            'total_sequences': sum(len(info['sequences']) for info in clusters.values())
        },
        'clusters': clusters,
        'tree_structure': tree_info,
        'weights': weights_info
    }
    
    # Write cluster analysis
    with open(output_dir / "cluster_analysis.txt", 'w') as f:
        f.write("WITCH Cluster Analysis Report\n")
        f.write("===========================\n\n")
        
        f.write("Summary:\n")
        f.write(f"- Total clusters: {report['summary']['total_clusters']}\n")
        f.write(f"- Total sequences: {report['summary']['total_sequences']}\n\n")
        
        f.write("Cluster Details:\n")
        for cluster_name, info in sorted(clusters.items()):
            f.write(f"\nCluster: {cluster_name}\n")
            f.write(f"- Number of sequences: {len(info['sequences'])}\n")
            f.write(f"- Has alignment: {info['has_alignment']}\n")
            f.write(f"- Has HMM: {info['has_hmm']}\n")
            f.write("- Sequences:\n")
            for seq in sorted(info['sequences']):
                f.write(f"  * {seq}\n")
        
        if tree_info:
            f.write("\nTree Structure:\n")
            f.write(tree_info['raw_tree'])
    
    # Write weights analysis to separate file
    if weights_info:
        write_weights_analysis(weights_info, output_dir / "weights_analysis.txt")
    
    # Save JSON for programmatic analysis
    with open(output_dir / "cluster_analysis.json", 'w') as f:
        json.dump(report, f, indent=2)

def main():
    if len(sys.argv) != 2:
        print("Usage: python analyze_witch_clusters.py <witch_output_dir>")
        sys.exit(1)
    
    analyze_witch_output(sys.argv[1])

if __name__ == "__main__":
    main() 