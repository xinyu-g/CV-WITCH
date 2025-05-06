#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Go up one level from src to get the project root
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Create output directory for evaluation results
mkdir -p "$PROJECT_ROOT/evaluation_results"

# Function to evaluate a dataset
evaluate_dataset() {
    local test_dir="$1"
    local ref_dir="$2"
    local output_dir="$3"
    local dataset_name="$4"
    local is_protein="$5"
    
    echo "Evaluating $dataset_name dataset..."
    echo "Test directory: $test_dir"
    echo "Reference directory: $ref_dir"
    echo "Output directory: $output_dir"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Run evaluation script
    if [ "$is_protein" = true ]; then
        python "$SCRIPT_DIR/evaluate_alignments.py" \
            --test-dir "$test_dir" \
            --ref-dir "$ref_dir" \
            --output-dir "$output_dir" \
            --protein
    else
        python "$SCRIPT_DIR/evaluate_alignments.py" \
            --test-dir "$test_dir" \
            --ref-dir "$ref_dir" \
            --output-dir "$output_dir"
    fi
    
    echo "----------------------------------------"
}

# Evaluate UPP results
echo "Evaluating UPP alignments..."
evaluate_dataset \
    "$PROJECT_ROOT/upp_output_benchmark/balibase" \
    "$PROJECT_ROOT/data/Outputs/balibase" \
    "$PROJECT_ROOT/evaluation_results/upp/balibase" \
    "BAliBASE" \
    true

evaluate_dataset \
    "$PROJECT_ROOT/upp_output_benchmark/1000M1" \
    "$PROJECT_ROOT/data/Outputs/1000M1" \
    "$PROJECT_ROOT/evaluation_results/upp/1000M1" \
    "1000M1" \
    false

evaluate_dataset \
    "$PROJECT_ROOT/upp_output_benchmark/1000M4" \
    "$PROJECT_ROOT/data/Outputs/1000M4" \
    "$PROJECT_ROOT/evaluation_results/upp/1000M4" \
    "1000M4" \
    false

# Evaluate WITCH results
echo "Evaluating WITCH alignments..."
evaluate_dataset \
    "$PROJECT_ROOT/witch_output_benchmark/balibase" \
    "$PROJECT_ROOT/data/Outputs/balibase" \
    "$PROJECT_ROOT/evaluation_results/witch/balibase" \
    "BAliBASE" \
    true

evaluate_dataset \
    "$PROJECT_ROOT/witch_output_benchmark/1000M1" \
    "$PROJECT_ROOT/data/Outputs/1000M1" \
    "$PROJECT_ROOT/evaluation_results/witch/1000M1" \
    "1000M1" \
    false

evaluate_dataset \
    "$PROJECT_ROOT/witch_output_benchmark/1000M4" \
    "$PROJECT_ROOT/data/Outputs/1000M4" \
    "$PROJECT_ROOT/evaluation_results/witch/1000M4" \
    "1000M4" \
    false

# Evaluate CV-enhanced WITCH results
echo "Evaluating CV-enhanced WITCH alignments..."
evaluate_dataset \
    "$PROJECT_ROOT/cv_witch_output/balibase" \
    "$PROJECT_ROOT/data/Outputs/balibase" \
    "$PROJECT_ROOT/evaluation_results/cv_witch/balibase" \
    "BAliBASE" \
    true

evaluate_dataset \
    "$PROJECT_ROOT/cv_witch_output/1000M1" \
    "$PROJECT_ROOT/data/Outputs/1000M1" \
    "$PROJECT_ROOT/evaluation_results/cv_witch/1000M1" \
    "1000M1" \
    false

evaluate_dataset \
    "$PROJECT_ROOT/cv_witch_output/1000M4" \
    "$PROJECT_ROOT/data/Outputs/1000M4" \
    "$PROJECT_ROOT/evaluation_results/cv_witch/1000M4" \
    "1000M4" \
    false

# Create a summary of all results
echo "Creating summary of all results..."
summary_file="$PROJECT_ROOT/evaluation_results/all_methods_summary.txt"

echo "Method\tDataset\tSPFN\tSPFP\tFN" > "$summary_file"

# Function to extract average metrics
get_avg_metrics() {
    local results_dir="$1"
    if [ -f "$results_dir/average_metrics.txt" ]; then
        spfn=$(grep "SPFN" "$results_dir/average_metrics.txt" | cut -d' ' -f3)
        spfp=$(grep "SPFP" "$results_dir/average_metrics.txt" | cut -d' ' -f3)
        fn=$(grep "FN" "$results_dir/average_metrics.txt" | cut -d' ' -f3)
        echo -e "${spfn}\t${spfp}\t${fn}"
    else
        echo -e "N/A\tN/A\tN/A"
    fi
}

# Add results for each method and dataset
for method in "upp" "witch" "cv_witch"; do
    for dataset in "balibase" "1000M1" "1000M4"; do
        results_dir="$PROJECT_ROOT/evaluation_results/${method}/${dataset}"
        metrics=$(get_avg_metrics "$results_dir")
        echo -e "${method}\t${dataset}\t${metrics}" >> "$summary_file"
    done
done

echo "Evaluation complete! Summary available at: $summary_file" 