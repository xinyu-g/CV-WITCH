#!/usr/bin/env python3

import os
import sys
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from sklearn.model_selection import KFold
import json
from collections import defaultdict
import shutil
from scipy.special import softmax

class CVEnhancedWitch:
    def __init__(self, witch_output_dir, k_folds=5, num_hmms=10, weighting_scheme='softmax'):
        """
        Initialize CV-enhanced WITCH processor
        
        Args:
            witch_output_dir: Directory containing WITCH output
            k_folds: Number of folds for cross-validation
            num_hmms: Number of top HMMs to use
            weighting_scheme: How to weight HMMs ('softmax', 'linear', or 'rank')
        """
        self.witch_dir = Path(witch_output_dir)
        self.k_folds = k_folds
        self.num_hmms = num_hmms
        self.weighting_scheme = weighting_scheme
        self.cv_scores = {}  # Store CV scores for each HMM
        
    def _compute_weights(self, scores):
        """
        Compute weights for HMMs based on CV scores using different schemes
        
        Args:
            scores: List of (hmm_id, score) tuples
        Returns:
            Dictionary of hmm_id to weight
        """
        if not scores:
            return {}
            
        hmm_ids, raw_scores = zip(*scores)
        raw_scores = np.array(raw_scores)
        
        if self.weighting_scheme == 'softmax':
            # Temperature parameter for softmax (adjust for sharper/smoother distribution)
            temp = 1.0
            weights = softmax(raw_scores / temp)
        elif self.weighting_scheme == 'linear':
            # Min-max normalization
            min_score = np.min(raw_scores)
            max_score = np.max(raw_scores)
            weights = (raw_scores - min_score) / (max_score - min_score + 1e-10)
        else:  # rank-based
            ranks = np.arange(len(raw_scores), 0, -1)
            weights = ranks / np.sum(ranks)
        
        return dict(zip(hmm_ids, weights))
    
    def _combine_alignments(self, alignments, weights):
        """
        Combine multiple alignments using weighted consensus
        
        Args:
            alignments: List of (alignment_file, weight) tuples
            weights: Dictionary of alignment weights
        Returns:
            Path to consensus alignment file
        """
        # Read all alignments
        alignment_data = []
        for aln_file, _ in alignments:
            seqs = self._read_fasta(aln_file)
            alignment_data.append(seqs)
        
        if not alignment_data:
            return None
            
        # Get all sequence IDs
        all_seq_ids = set()
        for aln in alignment_data:
            all_seq_ids.update(aln.keys())
        
        # Create consensus alignment
        consensus = {}
        for seq_id in all_seq_ids:
            # Get all versions of this sequence
            seq_versions = []
            seq_weights = []
            for (aln_file, hmm_id), aln in zip(alignments, alignment_data):
                if seq_id in aln:
                    seq_versions.append(aln[seq_id])
                    seq_weights.append(weights[hmm_id])
            
            if seq_versions:
                # Normalize weights for this sequence
                seq_weights = np.array(seq_weights)
                seq_weights = seq_weights / np.sum(seq_weights)
                
                # Use highest weighted version
                best_idx = np.argmax(seq_weights)
                consensus[seq_id] = seq_versions[best_idx]
        
        # Write consensus alignment
        consensus_file = Path(alignments[0][0]).parent / "consensus.afa"
        self._write_fasta(consensus, consensus_file)
        return consensus_file
    
    def process_cluster(self, cluster_dir):
        """
        Process a single cluster directory, performing CV on its sequences
        
        Args:
            cluster_dir: Path to cluster directory
        """
        cluster_id = cluster_dir.name
        input_fasta = cluster_dir / "hmmbuild.input.{}.fasta".format(cluster_id)
        
        if not input_fasta.exists():
            print(f"Warning: No input FASTA found for cluster {cluster_id}")
            return None
            
        # Read sequences
        sequences = self._read_fasta(input_fasta)
        if not sequences:
            return None
            
        # Perform k-fold CV
        kf = KFold(n_splits=self.k_folds, shuffle=True, random_state=42)
        fold_scores = []
        
        for fold_idx, (train_idx, val_idx) in enumerate(kf.split(list(sequences.items()))):
            # Create training and validation sets
            train_seqs = {k: sequences[k] for k in dict(list(sequences.items())[i] for i in train_idx)}
            val_seqs = {k: sequences[k] for k in dict(list(sequences.items())[i] for i in val_idx)}
            
            # Create temporary files
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_dir = Path(temp_dir)
                train_file = temp_dir / "train.fasta"
                val_file = temp_dir / "val.fasta"
                hmm_file = temp_dir / "model.hmm"
                
                # Write sequences to files
                self._write_fasta(train_seqs, train_file)
                self._write_fasta(val_seqs, val_file)
                
                # Build HMM from training set
                self._build_hmm(train_file, hmm_file)
                
                # Score validation sequences
                scores = self._score_sequences(hmm_file, val_file)
                if scores:
                    fold_scores.append(np.mean(list(scores.values())))
        
        # Calculate average CV score
        if fold_scores:
            avg_score = np.mean(fold_scores)
            self.cv_scores[cluster_id] = {
                'avg_score': avg_score,
                'fold_scores': fold_scores
            }
            return avg_score
        return None
    
    def _read_fasta(self, fasta_file):
        """Read sequences from FASTA file"""
        sequences = {}
        current_seq = []
        current_id = None
        
        with open(fasta_file) as f:
            for line in f:
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                        current_seq = []
                    current_id = line.strip()[1:]
                else:
                    current_seq.append(line.strip())
            if current_id:
                sequences[current_id] = ''.join(current_seq)
        return sequences
    
    def _write_fasta(self, sequences, output_file):
        """Write sequences to FASTA file"""
        with open(output_file, 'w') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n{seq}\n")
    
    def _build_hmm(self, input_file, output_file):
        """Build HMM from alignment"""
        cmd = ['hmmbuild', '--rna', str(output_file), str(input_file)]
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError:
            return False
    
    def _score_sequences(self, hmm_file, seq_file):
        """Score sequences against HMM"""
        cmd = ['hmmsearch', '--noali', '--tblout', '/dev/stdout', str(hmm_file), str(seq_file)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            scores = {}
            for line in result.stdout.split('\n'):
                if line and not line.startswith('#'):
                    fields = line.split()
                    seq_id = fields[0]
                    bit_score = float(fields[5])
                    scores[seq_id] = bit_score
            return scores
        except subprocess.CalledProcessError:
            return None
    
    def process_all_clusters(self):
        """Process all clusters in the WITCH output"""
        root_dir = self.witch_dir / "tree_decomp" / "root"
        if not root_dir.exists():
            raise FileNotFoundError(f"Root directory not found: {root_dir}")
        
        # Process each cluster
        for cluster_dir in root_dir.glob("A_*"):
            if cluster_dir.is_dir():
                print(f"Processing cluster {cluster_dir.name}...")
                self.process_cluster(cluster_dir)
        
        # Save CV scores
        output_file = self.witch_dir / "cv_scores.json"
        with open(output_file, 'w') as f:
            json.dump(self.cv_scores, f, indent=2)
    
    def align_with_cv_weights(self, input_fasta, output_file):
        """
        Align sequences using CV-weighted HMMs
        
        Args:
            input_fasta: Input sequences to align
            output_file: Where to save the final alignment
        """
        if not self.cv_scores:
            raise ValueError("No CV scores available. Run process_all_clusters first.")
        
        # Sort HMMs by CV score
        sorted_hmms = sorted(
            [(hmm_id, scores['avg_score']) for hmm_id, scores in self.cv_scores.items()],
            key=lambda x: x[1],
            reverse=True
        )
        
        # Select top HMMs
        top_hmms = sorted_hmms[:self.num_hmms]
        
        # Compute weights for selected HMMs
        weights = self._compute_weights(top_hmms)
        
        # Create temporary directory for intermediate files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            alignments = []
            
            # Align with each top HMM
            for hmm_id, _ in top_hmms:
                hmm_file = self.witch_dir / "tree_decomp" / "root" / hmm_id / f"hmmbuild.model.{hmm_id}"
                if hmm_file.exists():
                    aln_file = temp_dir / f"aln_{hmm_id}.afa"
                    
                    # Align sequences
                    cmd = ['hmmalign', '--outformat', 'afa', '-o', str(aln_file), str(hmm_file), input_fasta]
                    try:
                        subprocess.run(cmd, check=True)
                        alignments.append((aln_file, hmm_id))
                    except subprocess.CalledProcessError:
                        print(f"Warning: Failed to align with HMM {hmm_id}")
            
            if alignments:
                # Combine alignments using weights
                consensus_aln = self._combine_alignments(alignments, weights)
                if consensus_aln:
                    shutil.copy2(consensus_aln, output_file)
                    
                    # Save weights information
                    weights_file = Path(output_file).parent / "alignment_weights.json"
                    with open(weights_file, 'w') as f:
                        json.dump({
                            'weighting_scheme': self.weighting_scheme,
                            'weights': weights,
                            'cv_scores': {hmm_id: score for hmm_id, score in top_hmms}
                        }, f, indent=2)
                    return True
        return False

def main():
    if len(sys.argv) != 4:
        print("Usage: python cv_enhanced_witch.py <witch_output_dir> <input_fasta> <output_file>")
        sys.exit(1)
    
    witch_dir = sys.argv[1]
    input_fasta = sys.argv[2]
    output_file = sys.argv[3]
    
    cv_witch = CVEnhancedWitch(
        witch_dir,
        k_folds=5,
        num_hmms=10,
        weighting_scheme='softmax'
    )
    cv_witch.process_all_clusters()
    cv_witch.align_with_cv_weights(input_fasta, output_file)

if __name__ == "__main__":
    main() 