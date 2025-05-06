#!/usr/bin/env python3

import os
import sys
import subprocess
import random
import numpy as np
from pathlib import Path
from sklearn.model_selection import KFold
import shutil
import tempfile
import glob
import json

class WitchEnsembleCV:
    def __init__(self, input_fasta, k_folds=5, num_hmms=10, output_dir="witch_cv_output"):
        """
        Initialize the cross-validation ensemble builder
        
        Args:
            input_fasta: Path to input sequences
            k_folds: Number of folds for cross-validation
            num_hmms: Number of HMMs to use in WITCH
            output_dir: Directory for output files
        """
        self.input_fasta = input_fasta
        self.k_folds = k_folds
        self.num_hmms = num_hmms
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create directories for HMMs and alignments
        self.hmm_dir = self.output_dir / "hmms"
        self.hmm_dir.mkdir(exist_ok=True)
        self.aln_dir = self.output_dir / "alignments"
        self.aln_dir.mkdir(exist_ok=True)
        
    def split_sequences(self, fasta_file, fold_idx, temp_dir):
        """
        Split sequences into training and validation sets for a given fold
        
        Args:
            fasta_file: Input FASTA file
            fold_idx: Current fold index
            temp_dir: Directory for temporary files
        
        Returns:
            Paths to training and validation files
        """
        # Read sequences
        sequences = []
        with open(fasta_file) as f:
            current_seq = []
            current_header = None
            for line in f:
                if line.startswith('>'):
                    if current_header:
                        sequences.append((current_header, ''.join(current_seq)))
                        current_seq = []
                    current_header = line.strip()
                else:
                    current_seq.append(line.strip())
            if current_header:
                sequences.append((current_header, ''.join(current_seq)))
        
        # Create k-fold splits
        kf = KFold(n_splits=self.k_folds, shuffle=True, random_state=42)
        splits = list(kf.split(sequences))
        train_idx, val_idx = splits[fold_idx]
        
        # Create training and validation files
        train_file = os.path.join(temp_dir, f'train_fold_{fold_idx}.fasta')
        val_file = os.path.join(temp_dir, f'val_fold_{fold_idx}.fasta')
        
        with open(train_file, 'w') as f_train, open(val_file, 'w') as f_val:
            for idx in train_idx:
                header, seq = sequences[idx]
                f_train.write(f'{header}\n{seq}\n')
            for idx in val_idx:
                header, seq = sequences[idx]
                f_val.write(f'{header}\n{seq}\n')
                
        return train_file, val_file

    def build_overall_hmm(self, alignment_file, hmm_file):
        """
        Build a single HMM from a full alignment
        
        Args:
            alignment_file: Input alignment file
            hmm_file: Output HMM file
        """
        cmd = ['hmmbuild', '--rna', hmm_file, alignment_file]
        subprocess.run(cmd, check=True)

    def run_witch(self, input_file, output_dir):
        """
        Run WITCH on input sequences
        
        Args:
            input_file: Input FASTA file
            output_dir: Output directory for WITCH
        """
        # Create HMM directory within the WITCH output directory
        hmm_output_dir = Path(output_dir) / "hmms"
        hmm_output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            'witch-msa',
            '-i', input_file,
            '-d', output_dir,
            # '-p', str(hmm_output_dir),  # Directory to save HMMs
            '-k', str(self.num_hmms),
            '--keeptemp',  # Keep temporary files including HMMs
            '--save-weight', '1'  # Save HMM weights
        ]
        subprocess.run(cmd, check=True)

        # Build overall HMM from the final alignment
        if os.path.exists(os.path.join(output_dir, "aligned.fasta")):
            overall_hmm = os.path.join(hmm_output_dir, "overall.hmm")
            self.build_overall_hmm(os.path.join(output_dir, "aligned.fasta"), overall_hmm)

        return hmm_output_dir

    def score_sequences(self, hmm_file, sequences_file):
        """
        Score sequences against HMM using HMMER
        
        Args:
            hmm_file: HMM file
            sequences_file: Sequences to score
        
        Returns:
            Dictionary of sequence IDs to bit scores
        """
        cmd = ['hmmsearch', '--noali', '--tblout', '/dev/stdout', hmm_file, sequences_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        scores = {}
        for line in result.stdout.split('\n'):
            if line and not line.startswith('#'):
                fields = line.split()
                seq_id = fields[0]
                bit_score = float(fields[5])
                scores[seq_id] = bit_score
        return scores

    def align_with_hmm(self, hmm_file, sequences_file, output_file):
        """
        Align sequences using an HMM
        
        Args:
            hmm_file: HMM file to use
            sequences_file: Sequences to align
            output_file: Where to save the alignment
        """
        cmd = ['hmmalign', '--outformat', 'afa', '-o', output_file, hmm_file, sequences_file]
        subprocess.run(cmd, check=True)

    def run_cv(self):
        """
        Run the complete cross-validation pipeline
        """
        print("Starting cross-validation pipeline...")
        
        with tempfile.TemporaryDirectory() as temp_dir:
            fold_scores = {}  # Dictionary to store scores for each HMM across folds
            hmm_paths = {}    # Dictionary to store paths to HMMs
            
            # For each fold
            for fold in range(self.k_folds):
                print(f"\nProcessing fold {fold + 1}/{self.k_folds}")
                fold_dir = self.output_dir / f"fold_{fold}"
                fold_dir.mkdir(exist_ok=True)
                
                # Split data
                train_file, val_file = self.split_sequences(self.input_fasta, fold, temp_dir)
                
                # Run WITCH on training data
                witch_output = fold_dir / "witch_output"
                hmm_dir = self.run_witch(train_file, str(witch_output))
                
                # Get HMMs from WITCH decomposition
                hmm_files = glob.glob(str(hmm_dir / "*.hmm"))
                
                # Score validation sequences against each HMM
                for hmm_idx, hmm_file in enumerate(hmm_files):
                    hmm_name = f"hmm_{hmm_idx}"
                    if hmm_name not in fold_scores:
                        fold_scores[hmm_name] = []
                        hmm_paths[hmm_name] = []
                    
                    scores = self.score_sequences(hmm_file, val_file)
                    avg_score = np.mean(list(scores.values()))
                    fold_scores[hmm_name].append(avg_score)
                    hmm_paths[hmm_name].append(hmm_file)
                
                # Clean up temporary files
                os.remove(train_file)
                os.remove(val_file)
            
            # Find the best HMM based on average score across folds
            best_hmm, best_score = self.select_best_hmm(fold_scores, hmm_paths)
            
            # Copy the best HMM to the output directory
            best_hmm_final = self.hmm_dir / "best_hmm.hmm"
            shutil.copy2(best_hmm, best_hmm_final)
            
            # Align all sequences using the best HMM
            final_aln = self.aln_dir / "final_alignment.afa"
            self.align_with_hmm(best_hmm_final, self.input_fasta, final_aln)
            
            # Save results
            self.analyze_results(fold_scores, best_hmm, best_score)

    def select_best_hmm(self, fold_scores, hmm_paths):
        """
        Select the best HMM based on average score across folds
        """
        avg_scores = {hmm: np.mean(scores) for hmm, scores in fold_scores.items()}
        best_hmm_name = max(avg_scores, key=avg_scores.get)
        best_score = avg_scores[best_hmm_name]
        
        # Get the HMM file from the fold with the highest score for this HMM
        best_fold_idx = np.argmax(fold_scores[best_hmm_name])
        best_hmm_path = hmm_paths[best_hmm_name][best_fold_idx]
        
        return best_hmm_path, best_score

    def analyze_results(self, fold_scores, best_hmm, best_score):
        """
        Analyze cross-validation results and output summary
        """
        results_file = self.output_dir / "cv_results.txt"
        with open(results_file, 'w') as f:
            f.write("Cross-validation Results Summary\n")
            f.write("===============================\n\n")
            
            # Write detailed scores for each HMM and fold
            for hmm_name, scores in fold_scores.items():
                f.write(f"{hmm_name}:\n")
                for fold_idx, score in enumerate(scores):
                    f.write(f"  Fold {fold_idx + 1}: {score:.2f}\n")
                avg_score = np.mean(scores)
                f.write(f"  Average Score: {avg_score:.2f}\n\n")
            
            # Write best HMM information
            f.write("\nBest HMM Summary\n")
            f.write("================\n")
            f.write(f"Best HMM: {os.path.basename(best_hmm)}\n")
            f.write(f"Average Score: {best_score:.2f}\n")
            f.write(f"Final HMM location: {self.hmm_dir}/best_hmm.hmm\n")
            f.write(f"Final alignment: {self.aln_dir}/final_alignment.afa\n")
        
        # Also save detailed results in JSON format for potential further analysis
        with open(self.output_dir / "cv_results.json", 'w') as f:
            json.dump({
                'fold_scores': fold_scores,
                'best_hmm': str(best_hmm),
                'best_score': best_score
            }, f, indent=2)

def main():
    if len(sys.argv) not in [2, 3]:
        print("Usage: python cv_witch_ensemble.py <input_fasta> [output_dir]")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) == 3 else "witch_cv_output"
    
    cv_pipeline = WitchEnsembleCV(input_fasta, output_dir=output_dir)
    cv_pipeline.run_cv()

if __name__ == "__main__":
    main() 