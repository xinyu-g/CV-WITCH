#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Go up one level from src to get the project root
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Create output directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/upp_output"

# Function to process files with a specific pattern in a directory
process_directory() {
    local search_dir="$1"
    local file_pattern="$2"
    local output_prefix="$3"
    
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
            filename=$(basename "$file")
            filename_no_ext="${filename%.*}"  # Remove any extension
            output_dir="$PROJECT_ROOT/upp_output/$output_prefix/$dir_path/$filename_no_ext"
            
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
            
            # Run UPP with fixed FASTA file
            echo "Running UPP on fixed FASTA file..."
            run_upp.py -s "$fixed_fasta" -d "$output_dir"
            
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

# Process CRW-HF directory
echo "Processing CRW-HF directory..."
process_directory "$PROJECT_ROOT/data/CRW-HF" "*.unaln.fasta" "crw"

# Process 10k-HF directory
echo "Processing 10k-HF directory..."
process_directory "$PROJECT_ROOT/data/10k-HF" "unaligned-full.fas" "10k"

echo "Processing complete!" 