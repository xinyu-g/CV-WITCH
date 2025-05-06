#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Go up one level from src to get the project root
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Function to extract sequences not in training set
extract_holdout_seqs() {
    local full_fasta="$1"
    local train_fasta="$2"
    local output_fasta="$3"
    
    # Check if input files exist
    if [ ! -f "$full_fasta" ]; then
        echo "Error: Original FASTA file not found: $full_fasta"
        return 1
    fi
    
    if [ ! -f "$train_fasta" ]; then
        echo "Error: Training FASTA file not found: $train_fasta"
        return 1
    fi
    
    # Create a list of sequence IDs in the training set
    grep "^>" "$train_fasta" | sed 's/^>//' > "${output_fasta}.train_ids"
    
    # Extract sequences not in training set
    awk -v ids="${output_fasta}.train_ids" '
    BEGIN {
        while ((getline < ids) > 0) {
            exclude[$0] = 1
        }
        in_seq = 0
        keep = 0
    }
    /^>/ {
        id = substr($0, 2)
        in_seq = 1
        keep = !(id in exclude)
        if (keep) print
        next
    }
    in_seq && keep {
        print
    }' "$full_fasta" > "$output_fasta"
    
    # Clean up
    rm "${output_fasta}.train_ids"
}

# Process each fold in the CV results
for cv_dir in "$PROJECT_ROOT"/cv_ensemble_results/*/; do
    if [ -d "$cv_dir" ]; then
        dataset_name=$(basename "$cv_dir" _cv)
        echo "Processing $dataset_name..."
        
        # Check if all folds are complete
        if [ ! -d "${cv_dir}/fold_4" ]; then
            echo "Skipping $dataset_name - incomplete folds (missing fold_4)"
            echo "----------------------------------------"
            continue
        fi
        
        # Get the original FASTA file
        original_fasta="$PROJECT_ROOT/data/CRW/${dataset_name}.rawfasta"
        if [ ! -f "$original_fasta" ]; then
            echo "Cannot find original FASTA file: $original_fasta"
            continue
        fi
        
        # Process each fold
        for i in {0..4}; do
            fold_dir="${cv_dir}/fold_${i}"
            if [ -d "$fold_dir" ]; then
                echo "  Processing fold $i..."
                
                # Get training FASTA and WITCH output directory
                witch_output="${fold_dir}/witch_output"
                train_fasta="${witch_output}/aligned.fasta"
                
                if [ ! -d "$witch_output" ]; then
                    echo "    Cannot find WITCH output directory: $witch_output"
                    continue
                fi
                
                if [ ! -f "$train_fasta" ]; then
                    echo "    Cannot find training FASTA file: $train_fasta"
                    continue
                fi
                
                # Check if this fold has already been processed
                holdout_dir="${fold_dir}/holdout"
                if [ -d "${holdout_dir}/alignments" ]; then
                    echo "    Skipping fold $i - already processed"
                    echo "    ----------------------------------------"
                    continue
                fi
                
                # Extract held-out sequences
                mkdir -p "$holdout_dir"
                holdout_fasta="${holdout_dir}/holdout.fasta"
                
                echo "    Extracting held-out sequences..."
                if ! extract_holdout_seqs "$original_fasta" "$train_fasta" "$holdout_fasta"; then
                    echo "    Error extracting held-out sequences"
                    continue
                fi
                
                # Check if any sequences were extracted
                if [ ! -s "$holdout_fasta" ]; then
                    echo "    No held-out sequences found"
                    continue
                fi
                
                # Align held-out sequences using trained HMMs
                echo "    Aligning held-out sequences..."
                python "$SCRIPT_DIR/align_with_witch_hmms.py" \
                    "$witch_output" \
                    "$holdout_fasta" \
                    "${holdout_dir}/alignments"
                
                if [ $? -eq 0 ]; then
                    echo "    Successfully aligned held-out sequences"
                else
                    echo "    Error aligning held-out sequences"
                fi
                
                echo "    Done with fold $i"
                echo "----------------------------------------"
            fi
        done
    fi
done

echo "All processing complete!" 