#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Go up one level from src to get the project root
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Create output directory for CV-enhanced results
mkdir -p "$PROJECT_ROOT/cv_enhanced_witch_output"

# Function to process a single dataset
process_dataset() {
    local dataset_path="$1"
    local dataset_name=$(basename "$dataset_path")
    local output_dir="$PROJECT_ROOT/cv_enhanced_witch_output/$dataset_name"
    
    echo "----------------------------------------"
    echo "Processing dataset: $dataset_name"
    echo "Output directory: $output_dir"
    
    # Create dataset-specific output directory
    mkdir -p "$output_dir"
    
    # Get the WITCH output directory (now includes true_unalign)
    witch_output="$dataset_path/true_unalign"
    if [ ! -d "$witch_output/tree_decomp/root" ]; then
        echo "Error: Cannot find WITCH output directory structure in $witch_output"
        return 1
    fi
    
    # Get the unaligned sequences file from the original data
    # For balibase, it should be in data/Outputs/balibase/RV100_*/true_unalign.txt
    unaligned_file="$PROJECT_ROOT/data/Outputs/$(echo $dataset_path | sed 's|'$PROJECT_ROOT'/witch_output_benchmark/||')/true_unalign.txt"
    if [ ! -f "$unaligned_file" ]; then
        echo "Error: Cannot find unaligned sequences file: $unaligned_file"
        return 1
    fi
    
    # Run CV-enhanced WITCH
    echo "Running CV-enhanced WITCH..."
    python "$SCRIPT_DIR/cv_enhanced_witch.py" \
        "$witch_output" \
        "$unaligned_file" \
        "$output_dir/cv_enhanced_alignment.fasta"
    
    if [ $? -eq 0 ]; then
        echo "Successfully processed $dataset_name"
    else
        echo "Error processing $dataset_name"
    fi
    
    echo "----------------------------------------"
}

# Process balibase directory
echo "Processing balibase directory..."
for dataset in "$PROJECT_ROOT/witch_output_benchmark/balibase"/RV100_*; do
    if [ -d "$dataset" ]; then
        process_dataset "$dataset"
    fi
done

# Process 1000M1 directory (R0-R9)
echo "Processing 1000M1 directory..."
for i in {0..9}; do
    for dataset in "$PROJECT_ROOT/witch_output_benchmark/1000M1/R$i"/*; do
        if [ -d "$dataset" ]; then
            process_dataset "$dataset"
        fi
    done
done

# Process 1000M4 directory (R0-R9)
echo "Processing 1000M4 directory..."
for i in {0..9}; do
    for dataset in "$PROJECT_ROOT/witch_output_benchmark/1000M4/R$i"/*; do
        if [ -d "$dataset" ]; then
            process_dataset "$dataset"
        fi
    done
done

echo "All processing complete!"

# Create a summary of results
echo "Creating summary of results..."
summary_file="$PROJECT_ROOT/cv_enhanced_witch_output/summary.txt"
echo "CV-Enhanced WITCH Results Summary - $(date)" > "$summary_file"
echo "----------------------------------------" >> "$summary_file"

# Group results by dataset type
for dataset_type in "balibase" "1000M1" "1000M4"; do
    echo "Results for $dataset_type:" >> "$summary_file"
    echo "=========================" >> "$summary_file"
    
    find "$PROJECT_ROOT/cv_enhanced_witch_output" -path "*/$dataset_type/*" -name "alignment_weights.json" | while read -r weights_file; do
        dataset_name=$(basename $(dirname "$weights_file"))
        echo "Dataset: $dataset_name" >> "$summary_file"
        echo "Weights and CV Scores:" >> "$summary_file"
        python -m json.tool "$weights_file" >> "$summary_file"
        echo "----------------------------------------" >> "$summary_file"
    done
done

echo "Summary saved to: $summary_file"

# Create a simple comparison script
cat > "$PROJECT_ROOT/cv_enhanced_witch_output/compare_alignments.sh" << 'EOF'
#!/bin/bash
# Compare original WITCH alignments with CV-enhanced alignments

PROJECT_ROOT="$(dirname "$(dirname "$0")")"

for cv_aln in $(find "$PROJECT_ROOT/cv_enhanced_witch_output" -name "cv_enhanced_alignment.fasta"); do
    dataset_path=$(dirname "$cv_aln")
    dataset_name=$(basename "$dataset_path")
    original_aln="$PROJECT_ROOT/witch_output_benchmark/$(echo $dataset_path | sed 's|'$PROJECT_ROOT'/cv_enhanced_witch_output/||')/true_unalign/aligned.fasta"
    
    if [ -f "$original_aln" ]; then
        echo "Comparing alignments for $dataset_name..."
        # Add your comparison metrics here
        # For example: qscore, sp-score, etc.
    fi
done
EOF

chmod +x "$PROJECT_ROOT/cv_enhanced_witch_output/compare_alignments.sh" 