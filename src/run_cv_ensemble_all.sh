#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Go up one level from src to get the project root
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Create output directory for results if it doesn't exist
mkdir -p "$PROJECT_ROOT/cv_ensemble_results"

# Function to process a single dataset
process_dataset() {
    local dataset="$1"
    local dataset_name=$(basename "$dataset")
    local output_dir="$PROJECT_ROOT/cv_ensemble_results/$dataset_name"
    
    echo "----------------------------------------"
    echo "Processing dataset: $dataset_name"
    echo "Output directory: $output_dir"
    
    # Create dataset-specific output directory
    mkdir -p "$output_dir"
    
    # Check if tree_decomp/root exists for this dataset
    if [ ! -d "$dataset/tree_decomp/root" ]; then
        echo "Error: tree_decomp/root not found for $dataset_name"
        return 1
    fi
    
    # Run the Python script for this dataset
    echo "Running cv_hmm_ensemble.py on $dataset_name..."
    python3 "$SCRIPT_DIR/cv_hmm_ensemble.py" "$dataset/tree_decomp/root"
    
    # Check if the run was successful
    if [ $? -eq 0 ]; then
        # Move the results to the dataset-specific output directory
        if [ -f "cv_results.txt" ]; then
            mv cv_results.txt "$output_dir/"
            echo "Results saved to $output_dir/cv_results.txt"
        fi
        if [ -f "cv_hmm_ensemble.log" ]; then
            mv cv_hmm_ensemble.log "$output_dir/"
            echo "Log saved to $output_dir/cv_hmm_ensemble.log"
        fi
        echo "Successfully processed $dataset_name"
    else
        echo "Error processing $dataset_name"
    fi
    
    echo "----------------------------------------"
}

# Main execution
echo "Starting CV HMM ensemble processing for all datasets..."

# Process each dataset in the witch_output directory
for dataset in "$PROJECT_ROOT/witch_output/16S."*; do
    if [ -d "$dataset" ]; then
        process_dataset "$dataset"
    fi
done

echo "Processing complete!"

# Create a summary of all results
echo "Creating summary of all results..."
summary_file="$PROJECT_ROOT/cv_ensemble_results/all_datasets_summary.txt"
echo "Dataset Summary - $(date)" > "$summary_file"
echo "----------------------------------------" >> "$summary_file"

for result_dir in "$PROJECT_ROOT/cv_ensemble_results/16S."*; do
    if [ -d "$result_dir" ]; then
        dataset_name=$(basename "$result_dir")
        echo "Results for $dataset_name:" >> "$summary_file"
        if [ -f "$result_dir/cv_results.txt" ]; then
            echo "Top performing clusters:" >> "$summary_file"
            head -n 5 "$result_dir/cv_results.txt" >> "$summary_file"
        else
            echo "No results found" >> "$summary_file"
        fi
        echo "----------------------------------------" >> "$summary_file"
    fi
done

echo "Summary saved to: $summary_file" 