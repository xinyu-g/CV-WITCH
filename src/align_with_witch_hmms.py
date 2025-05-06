#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path
import tempfile
from collections import defaultdict
import numpy as np
import re

class WitchHMMAligner:
    def __init__(self, witch_output_dir, num_hmms=10):
        """
        Initialize aligner with WITCH output directory containing HMMs and weights
        
        Args:
            witch_output_dir: Directory containing WITCH output (HMMs and weights)
            num_hmms: Number of top HMMs to use for each sequence
        """
        self.witch_dir = Path(witch_output_dir)
        self.num_hmms = num_hmms
        self.hmm_files = self._find_hmms()
        # Make weights optional
        try:
            self.weights = self._load_weights()
        except (FileNotFoundError, ValueError) as e:
            print(f"Warning: Could not load weights file ({e}). Proceeding without weights.")
            self.weights = None
        
    def _load_weights(self):
        """Load and parse weights file"""
        weights_file = self.witch_dir / "weights.txt"
        if not weights_file.exists():
            return None
            
        weights = defaultdict(dict)  # seq_id -> hmm_id -> weight
        with open(weights_file) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    try:
                        # Try different parsing approaches
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            # Try to handle different formats
                            seq_id = parts[0].strip("(',)")
                            hmm_id = parts[1].strip("(',)")
                            # Find the first float-like string
                            weight = None
                            for part in parts[2:]:
                                try:
                                    weight = float(part.strip("(',)"))
                                    break
                                except ValueError:
                                    continue
                            if weight is not None:
                                weights[seq_id][hmm_id] = weight
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Skipping malformed line in weights file: {line.strip()}")
                        continue
        return weights if weights else None
    
    def _find_hmms(self):
        """Find all HMM files in the WITCH output"""
        hmm_files = {}
        
        # Check tree_decomp/root directory for HMM files
        root_dir = self.witch_dir / "tree_decomp" / "root"
        if not root_dir.exists():
            raise FileNotFoundError(f"Root directory not found: {root_dir}")
            
        # Look for cluster directories (A_*) and their HMM model files
        for cluster_dir in root_dir.glob("A_*"):
            if cluster_dir.is_dir():
                cluster_id = cluster_dir.name
                hmm_file = cluster_dir / f"hmmbuild.model.{cluster_id}"
                if hmm_file.exists():
                    hmm_files[cluster_id] = hmm_file
                else:
                    print(f"Warning: No HMM file found for cluster {cluster_id}")
        
        if not hmm_files:
            raise FileNotFoundError(f"No HMM files found in {root_dir}")
            
        print(f"Found {len(hmm_files)} HMM files")
        return hmm_files
    
    def score_sequence(self, seq_file, hmm_file):
        """Score a sequence against an HMM using hmmsearch"""
        cmd = ['hmmsearch', '--noali', '--tblout', '/dev/stdout', str(hmm_file), seq_file]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        scores = {}
        for line in result.stdout.split('\n'):
            if line and not line.startswith('#'):
                try:
                    fields = line.split()
                    seq_id = fields[0]
                    bit_score = float(fields[5])
                    scores[seq_id] = bit_score
                except (IndexError, ValueError):
                    continue
        return scores
    
    def align_with_hmm(self, seq_file, hmm_file, output_file):
        """Align sequences using hmmalign"""
        cmd = ['hmmalign', '--outformat', 'afa', '-o', output_file, str(hmm_file), seq_file]
        subprocess.run(cmd, check=True)
    
    def align_sequences(self, input_fasta, output_dir):
        """
        Align sequences using the trained HMMs and weights
        
        Args:
            input_fasta: FASTA file containing sequences to align
            output_dir: Directory to save results
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create temporary directory for intermediate files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            
            # Score sequences against all HMMs
            print("Scoring sequences against HMMs...")
            seq_scores = defaultdict(dict)  # seq_id -> hmm_id -> score
            for hmm_id, hmm_file in self.hmm_files.items():
                print(f"  Scoring against {hmm_id}...")
                scores = self.score_sequence(input_fasta, hmm_file)
                for seq_id, score in scores.items():
                    seq_scores[seq_id][hmm_id] = score
            
            # For each sequence, select top HMMs and align
            print("Aligning sequences...")
            for seq_id in seq_scores:
                # Sort HMMs by score for this sequence
                hmm_scores = [(hmm_id, score) for hmm_id, score in seq_scores[seq_id].items()]
                top_hmms = sorted(hmm_scores, key=lambda x: x[1], reverse=True)[:self.num_hmms]
                
                # Align with each top HMM
                alignments = []
                for hmm_id, score in top_hmms:
                    aln_file = temp_dir / f"{seq_id}_{hmm_id}.afa"
                    self.align_with_hmm(input_fasta, self.hmm_files[hmm_id], aln_file)
                    alignments.append((aln_file, score))
                
                # Combine alignments (for now, just use the best one)
                best_aln = max(alignments, key=lambda x: x[1])[0]
                output_file = output_dir / f"{seq_id}_aligned.afa"
                with open(best_aln) as f_in, open(output_file, 'w') as f_out:
                    f_out.write(f_in.read())
            
            # Write alignment scores
            scores_file = output_dir / "alignment_scores.txt"
            with open(scores_file, 'w') as f:
                f.write("Sequence\tHMM\tScore\n")
                for seq_id in seq_scores:
                    for hmm_id, score in seq_scores[seq_id].items():
                        f.write(f"{seq_id}\t{hmm_id}\t{score}\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python align_with_witch_hmms.py <witch_output_dir> <input_fasta> <output_dir>")
        sys.exit(1)
    
    witch_dir = sys.argv[1]
    input_fasta = sys.argv[2]
    output_dir = sys.argv[3]
    
    aligner = WitchHMMAligner(witch_dir)
    aligner.align_sequences(input_fasta, output_dir)

if __name__ == "__main__":
    main() 