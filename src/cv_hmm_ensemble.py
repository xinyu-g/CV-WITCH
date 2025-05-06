#!/usr/bin/env python3

import os
import sys
import glob
import numpy as np
from sklearn.model_selection import KFold
import subprocess
import logging
from pathlib import Path
import shutil
import tempfile
import re

def setup_logging():
    """Configure logging for the script"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('cv_hmm_ensemble.log')
        ]
    )

def clean_sequence(seq, keep_gaps=True):
    """Clean a sequence, optionally preserving alignment gaps.
    
    Args:
        seq (str): Input sequence
        keep_gaps (bool): Whether to preserve alignment gaps
        
    Returns:
        str: Cleaned sequence
    """
    if keep_gaps:
        # Keep alignment gaps (-) and valid nucleotides
        valid_chars = set('ACGTUacgtu-')
        cleaned = ''.join(c if c in valid_chars else '-' for c in seq)
    else:
        # Remove gaps and keep only valid nucleotides
        valid_chars = set('ACGTUacgtu')
        cleaned = ''.join(c for c in seq if c in valid_chars)
    return cleaned.upper()

def read_fasta(filename):
    """Read sequences from a FASTA file.
    
    Args:
        filename (str): Path to the FASTA file
        
    Returns:
        list: List of (header, sequence) tuples with cleaned sequences
    """
    sequences = []
    current_header = None
    current_seq = []
    expected_length = None
    
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    seq = ''.join(current_seq)
                    cleaned_seq = clean_sequence(seq)
                    # Set expected length from first sequence
                    if expected_length is None:
                        expected_length = len(cleaned_seq)
                        logging.info(f"Expected sequence length: {expected_length}")
                    # Verify sequence length
                    if len(cleaned_seq) != expected_length:
                        logging.error(f"Sequence length mismatch for {current_header}: "
                                   f"got {len(cleaned_seq)}, expected {expected_length}")
                        continue
                    sequences.append((current_header, cleaned_seq))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
                
    if current_header:
        seq = ''.join(current_seq)
        cleaned_seq = clean_sequence(seq)
        if expected_length is None or len(cleaned_seq) == expected_length:
            sequences.append((current_header, cleaned_seq))
        else:
            logging.error(f"Sequence length mismatch for {current_header}: "
                       f"got {len(cleaned_seq)}, expected {expected_length}")
    
    # Log sequence cleaning statistics
    if sequences:
        lengths = [len(seq) for _, seq in sequences]
        if len(set(lengths)) > 1:
            logging.error(f"Found sequences of different lengths: {set(lengths)}")
            return []
        logging.info(f"Read {len(sequences)} valid sequences from {filename}, "
                    f"all with length {lengths[0]}")
    else:
        logging.warning(f"No valid sequences found in {filename} after cleaning")
    
    return sequences

def write_fasta(sequences, filename, keep_gaps=True):
    """Write sequences to a FASTA file.
    
    Args:
        sequences (list): List of (header, sequence) tuples
        filename (str): Output file path
        keep_gaps (bool): Whether to preserve alignment gaps
    """
    try:
        # Verify all sequences have the same length if keeping gaps
        if keep_gaps and sequences:
            lengths = [len(seq) for _, seq in sequences]
            if len(set(lengths)) > 1:
                logging.error(f"Cannot write sequences of different lengths: {set(lengths)}")
                return False
        
        valid_sequences = [(header, seq) for header, seq in sequences if seq]
        if len(valid_sequences) != len(sequences):
            logging.warning(f"Skipping {len(sequences) - len(valid_sequences)} empty sequences")
        
        if not valid_sequences:
            logging.error(f"No valid sequences to write to {filename}")
            return False
            
        with open(filename, 'w') as f:
            for header, seq in valid_sequences:
                # Clean sequence based on whether we're keeping gaps
                cleaned_seq = clean_sequence(seq, keep_gaps=keep_gaps)
                if not cleaned_seq:
                    logging.warning(f"Skipping empty sequence after cleaning for header: {header}")
                    continue
                    
                # Write in chunks of 60 characters for better readability
                f.write(f">{header}\n")
                for i in range(0, len(cleaned_seq), 60):
                    f.write(f"{cleaned_seq[i:i+60]}\n")
                f.flush()
        
        # Verify the file was written
        if not os.path.exists(filename):
            logging.error(f"Failed to create FASTA file: {filename}")
            return False
        if os.path.getsize(filename) == 0:
            logging.error(f"FASTA file is empty: {filename}")
            return False
            
        logging.info(f"Successfully wrote {len(valid_sequences)} sequences to {filename}")
        return True
    except Exception as e:
        logging.error(f"Error writing FASTA file {filename}: {str(e)}")
        return False

def train_hmm(fasta_file, output_model):
    """Train an HMM using HMMER's hmmbuild.
    
    Args:
        fasta_file (str): Input FASTA file path
        output_model (str): Output HMM model file path
        
    Returns:
        bool: True if training succeeded, False otherwise
    """
    try:
        # Add --cpu 1 to avoid threading issues
        # Add --fast for faster building during cross-validation
        result = subprocess.run(
            ['hmmbuild', '--cpu', '1', '--fast', output_model, fasta_file],
            check=True, capture_output=True, text=True
        )
        
        # Verify the model file was created and has content
        if not os.path.exists(output_model) or os.path.getsize(output_model) == 0:
            logging.error(f"HMM model file was not created or is empty: {output_model}")
            return False
            
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error training HMM: {e}")
        logging.error(f"STDERR: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Unexpected error in train_hmm: {str(e)}")
        return False

def evaluate_hmm(model_file, test_fasta, output_file):
    """Evaluate an HMM using HMMER's hmmsearch.
    
    Args:
        model_file (str): Path to the HMM model
        test_fasta (str): Path to the test sequences
        output_file (str): Path to save the results
        
    Returns:
        float: Average bit score across all test sequences
    """
    try:
        # Add --max to turn off heuristic filters and --noali to reduce output size
        # Add -E 1000 to ensure we get scores even for weak matches
        # Add --cpu 1 to avoid any threading issues
        result = subprocess.run(
            ['hmmsearch', '--max', '--noali', '-E', '1000', '--cpu', '1',
             '--tblout', output_file, model_file, test_fasta],
            check=True, capture_output=True, text=True
        )
        
        # Check if the output file exists and has content
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            logging.warning(f"No hits found in hmmsearch for {test_fasta}")
            return 0.0
        
        # Parse bit scores from the output file
        scores = []
        with open(output_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) >= 6:  # Ensure we have enough fields
                    try:
                        scores.append(float(fields[5]))  # Bit score is in 6th column
                    except (ValueError, IndexError) as e:
                        logging.warning(f"Could not parse score from line: {line.strip()}")
                        continue
                
        return np.mean(scores) if scores else 0.0
    except subprocess.CalledProcessError as e:
        logging.error(f"Error evaluating HMM: {e}")
        logging.error(f"STDERR: {e.stderr}")
        return 0.0
    except Exception as e:
        logging.error(f"Unexpected error in evaluate_hmm: {str(e)}")
        return 0.0

def process_cluster(cluster_path, k=5, min_seq_length=10):
    """Process a single cluster with k-fold cross-validation.
    
    Args:
        cluster_path (str): Path to the cluster directory
        k (int): Number of folds for cross-validation
        min_seq_length (int): Minimum sequence length to consider
        
    Returns:
        tuple: (best_model_path, avg_performance) or (None, 0) on failure
    """
    input_fasta = os.path.join(cluster_path, 'hmmbuild.input.*.fasta')
    input_files = glob.glob(input_fasta)
    
    if not input_files:
        logging.error(f"No input FASTA found in {cluster_path}")
        return None, 0
    
    # Read sequences from the first matching file
    input_file = input_files[0]
    logging.info(f"Reading sequences from {input_file}")
    sequences = read_fasta(input_file)
    
    if not sequences:
        logging.error(f"No sequences found in {input_file}")
        return None, 0
    
    # Filter sequences by length
    valid_sequences = [(header, seq) for header, seq in sequences if len(seq) >= min_seq_length]
    if len(valid_sequences) != len(sequences):
        logging.warning(f"Filtered out {len(sequences) - len(valid_sequences)} sequences shorter than {min_seq_length} bp")
        sequences = valid_sequences
    
    if len(sequences) < k:
        logging.error(f"Not enough sequences ({len(sequences)}) for {k}-fold cross-validation")
        return None, 0
    
    logging.info(f"Found {len(sequences)} valid sequences in {input_file}")
    
    # Create a persistent temporary directory
    cluster_name = os.path.basename(cluster_path)
    temp_dir = os.path.join(cluster_path, 'cv_temp')
    try:
        # Create temp directory if it doesn't exist
        os.makedirs(temp_dir, exist_ok=True)
        logging.info(f"Created temporary directory: {temp_dir}")
        
        kf = KFold(n_splits=k, shuffle=True, random_state=42)
        fold_models = []
        fold_scores = []
        
        # Perform k-fold cross-validation
        for fold_idx, (train_idx, test_idx) in enumerate(kf.split(sequences)):
            train_seqs = [sequences[i] for i in train_idx]
            test_seqs = [sequences[i] for i in test_idx]
            
            if not train_seqs or not test_seqs:
                logging.error(f"Empty split for fold {fold_idx + 1}")
                continue
            
            # Create temporary files for this fold
            train_file = os.path.join(temp_dir, f'train_fold_{fold_idx}.fasta')
            test_file = os.path.join(temp_dir, f'test_fold_{fold_idx}.fasta')
            model_file = os.path.join(temp_dir, f'model_fold_{fold_idx}.hmm')
            eval_file = os.path.join(temp_dir, f'eval_fold_{fold_idx}.txt')
            
            # Write sequences to files and verify they exist
            logging.info(f"Writing {len(train_seqs)} training sequences to {train_file}")
            if not write_fasta(train_seqs, train_file, keep_gaps=True):  # Keep gaps for hmmbuild
                logging.error(f"Failed to write training sequences to {train_file}")
                continue
                
            logging.info(f"Writing {len(test_seqs)} test sequences to {test_file}")
            if not write_fasta(test_seqs, test_file, keep_gaps=False):  # Remove gaps for hmmsearch
                logging.error(f"Failed to write test sequences to {test_file}")
                continue
            
            # Train and evaluate HMM
            if train_hmm(train_file, model_file):
                # Verify model file exists before evaluation
                if not os.path.exists(model_file) or os.path.getsize(model_file) == 0:
                    logging.error(f"HMM model file not found or empty: {model_file}")
                    continue
                    
                score = evaluate_hmm(model_file, test_file, eval_file)
                fold_models.append(model_file)
                fold_scores.append(score)
                logging.info(f"Fold {fold_idx + 1} score: {score}")
            else:
                logging.error(f"Failed to train HMM for fold {fold_idx + 1}")
        
        if not fold_scores:
            logging.error(f"No successful folds for cluster {os.path.basename(cluster_path)}")
            return None, 0
        
        # Select the best model based on average performance
        best_idx = np.argmax(fold_scores)
        best_score = fold_scores[best_idx]
        
        # Copy the best model to the cluster directory
        best_model = fold_models[best_idx]
        output_model = os.path.join(cluster_path, 'best_cv_model.hmm')
        try:
            shutil.copy2(best_model, output_model)
            if not os.path.exists(output_model) or os.path.getsize(output_model) == 0:
                logging.error(f"Failed to copy best model to {output_model}")
                return None, 0
            logging.info(f"Successfully saved best model (score: {best_score}) to {output_model}")
        except Exception as e:
            logging.error(f"Error copying best model: {str(e)}")
            return None, 0
        
        # Clean up temporary files if successful
        try:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                logging.info(f"Cleaned up temporary directory: {temp_dir}")
        except Exception as e:
            logging.warning(f"Failed to clean up temporary directory {temp_dir}: {str(e)}")
            
        return output_model, best_score
            
    except Exception as e:
        logging.error(f"Error in process_cluster: {str(e)}")
        # Try to clean up on error
        try:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
        except:
            pass
        return None, 0

def main():
    """Main function to process all clusters."""
    setup_logging()
    
    if len(sys.argv) != 2:
        print("Usage: python cv_hmm_ensemble.py <path_to_tree_decomp_root>")
        sys.exit(1)
    
    root_dir = sys.argv[1]
    cluster_pattern = os.path.join(root_dir, 'A_0_*')
    clusters = glob.glob(cluster_pattern)
    
    results = []
    for cluster in sorted(clusters):
        cluster_name = os.path.basename(cluster)
        logging.info(f"Processing cluster {cluster_name}")
        
        model_path, performance = process_cluster(cluster)
        if model_path:
            results.append((cluster_name, performance))
            logging.info(f"Cluster {cluster_name} - Average performance: {performance}")
        else:
            logging.error(f"Failed to process cluster {cluster_name}")
    
    # Write summary of results
    with open('cv_results.txt', 'w') as f:
        for cluster_name, performance in sorted(results, key=lambda x: x[1], reverse=True):
            f.write(f"{cluster_name}\t{performance}\n")

if __name__ == "__main__":
    main() 