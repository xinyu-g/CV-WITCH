#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Go up one level from src to get the project root
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Create output directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/upp_output_benchmark"

# Function to process files with a specific pattern in a directory
process_directory() {
    local search_dir="$1"
    local file_pattern="$2"
    local output_prefix="$3"
    local molecule_type="$4"  # Added molecule type parameter
    local dataset_type="$5"   # Added dataset type parameter
    
    # Create a sorted list of files by size (smallest to largest)
    local dir_files=()
    while IFS= read -r line; do
        dir_files+=("$line")
    done < <(find "$search_dir" -name "$file_pattern" -type f -exec ls -s {} \; | sort -n | awk '{print $2}')
    
    echo "Found ${#dir_files[@]} files matching $file_pattern in $search_dir"
    echo "Files in order of processing (smallest to largest):"
    for file in "${dir_files[@]}"; do
        echo "$(ls -sh "$file")"
    done
    echo "----------------------------------------"
    
    # Process each file in the sorted list
    for file in "${dir_files[@]}"; do
        if [ -f "$file" ]; then
            # Get the relative path structure
            rel_path=${file#"$search_dir/"}
            dir_path=$(dirname "$rel_path")
            filename=$(basename "$file" .txt)  # Remove .txt extension
            output_dir="$PROJECT_ROOT/upp_output_benchmark/$output_prefix/$dir_path/$filename"
            
            # Check if the directory exists and has completed alignment
            if [ -f "$output_dir/output_alignment.fasta" ]; then
                echo "Skipping $filename - output already exists and complete"
                echo "----------------------------------------"
                continue
            else
                # If directory exists but is incomplete, remove it
                if [ -d "$output_dir" ]; then
                    echo "Found incomplete output directory for $filename - removing and reprocessing"
                    rm -rf "$output_dir"
                fi
            fi
            
            echo "Processing contents of $file..."
            echo "File size: $(ls -sh "$file" | awk '{print $1}')"
            echo "Output will be saved to: $output_dir"
            
            # Create output directory (including parent directories)
            mkdir -p "$output_dir"
            
            # Fix FASTA formatting
            fixed_fasta="$output_dir/fixed.fasta"
            echo "Fixing FASTA format..."
            python "$SCRIPT_DIR/fix_fasta_format.py" "$file" "$fixed_fasta"
            
            if [ $? -ne 0 ]; then
                echo "Error fixing FASTA format for $filename"
                echo "----------------------------------------"
                continue
            fi
            
            # Set backbone size based on dataset type
            local backbone_size="100"  # Default for balibase
            if [[ "$dataset_type" == "1000M"* ]]; then
                # For 1000M datasets, use 20% of sequences as backbone
                total_seqs=$(grep -c "^>" "$fixed_fasta")
                backbone_size=$((total_seqs / 5))  # 20% of total sequences
                echo "Setting backbone size to $backbone_size (20% of $total_seqs sequences)"
            fi
            
            # Run UPP with fixed FASTA file and appropriate molecule type
            echo "Running UPP on fixed FASTA file (molecule type: $molecule_type, backbone size: $backbone_size)..."
            run_upp.py -s "$fixed_fasta" -m "$molecule_type" -d "$output_dir" -b "$backbone_size"
            
            # Check if the run was successful and output_alignment.fasta exists
            if [ $? -eq 0 ] && [ -f "$output_dir/output_alignment.fasta" ]; then
                echo "Successfully completed UPP for $filename"
                # Clean up fixed FASTA file
                rm "$fixed_fasta"
            else
                echo "Error processing $filename or output file not found"
            fi
            
            echo "----------------------------------------"
        fi
    done
}

# Process balibase directory (protein sequences)
echo "Processing balibase directory..."
process_directory "$PROJECT_ROOT/data/Outputs/balibase" "true_unalign.txt" "balibase" "amino" "balibase"

# Process 1000M1 directory (R0-R9) (RNA sequences)
echo "Processing 1000M1 directory..."
for i in {0..9}; do
    process_directory "$PROJECT_ROOT/data/Outputs/1000M1/R$i" "true_unalign.txt" "1000M1/R$i" "dna" "1000M1"
done

# Process 1000M4 directory (R0-R9) (RNA sequences)
echo "Processing 1000M4 directory..."
for i in {0..9}; do
    process_directory "$PROJECT_ROOT/data/Outputs/1000M4/R$i" "true_unalign.txt" "1000M4/R$i" "dna" "1000M4"
done

echo "Processing complete!" 