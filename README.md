# Benchmark Analysis: UPP vs WITCH vs CV-WITCH

This repository contains scripts and tools for running and comparing three multiple sequence alignment methods: UPP, WITCH, and CV-enhanced WITCH. Below are detailed instructions for running each method and generating benchmark comparisons.

## Directory Structure
```
.
├── data (not included here)/
│   ├── CRW-HF/
│   │   ├── 16S.3/
│   │   ├── 16S.B.ALL/
│   │   └── 16S.T/
│   ├── balibase/
│   ├── 1000M1/
│   └── 1000M4/
├── src/
│   ├── run_upp_balibase_1000M.sh
│   ├── run_witch_balibase_1000M.sh
│   ├── run_cv_enhanced_witch.sh
│   └── plot_benchmark_results.py
└── plots/
    ├── crw_benchmark.png
    ├── balibase_benchmark.png
    ├── 1000M1_benchmark.png
    └── 1000M4_benchmark.png
```

## Prerequisites

1. Install required software:
```bash
# UPP
See https://github.com/smirarab/sepp/blob/master/README.UPP.md

# WITCH dependencies
See https://github.com/c5shen/WITCH
```

## Running UPP

1. For CRW-HF datasets (see src/run_upp.sh):
```bash
# Process each dataset
./src/run_upp.sh
```

2. For BAliBASE and 1000M datasets:
```bash
# Run the UPP processing script
./src/run_upp_balibase_1000M.sh
```

The UPP results will be saved in `upp_output_benchmark/` directory.

## Running WITCH

1. For CRW-HF datasets:
```bash
# Process each dataset
./src/run_witch.sh
```

2. For BAliBASE and 1000M datasets:
```bash
# Run the WITCH processing script
./src/run_witch_balibase_1000M.sh
```

The WITCH results will be saved in `witch_output_benchmark/` directory.

## Running CV-enhanced WITCH

1. For all datasets:
```bash
# Run the CV-enhanced WITCH script
./src/run_cv_enhanced_witch.sh
```

This script will:
- Process each dataset using k-fold cross-validation (k=5)
- Generate HMM weights based on CV performance
- Create final alignments using weighted HMM consensus
- Save results in `cv_witch_output/` directory

## Generating Benchmark Plots

After running all three methods, generate comparison plots:
```bash
python src/plot_benchmark_results.py
```

This will create four plots in the `plots/` directory:
- `crw_benchmark.png`: Comparison on CRW-HF datasets
- `balibase_benchmark.png`: Comparison on BAliBASE datasets
- `1000M1_benchmark.png`: Comparison on 1000M1 datasets
- `1000M4_benchmark.png`: Comparison on 1000M4 datasets

## Key Parameters

### UPP Parameters
- Backbone size: 10% of sequences for  BAliBASE and 1000M datasets
- Molecule type: amino for BAliBASE, dna for 1000M datasets

### WITCH Parameters
- Number of HMMs: 10
- Decomposition size: 10% of sequences for  BAliBASE and 1000M datasets
- Save weights: enabled

### CV-WITCH Parameters
- Number of folds: 5
- Number of HMMs: 10
- Weighting scheme: Linear min-max
- Cross-validation scoring: bit score based

## Output Files

Each method generates the following outputs:

### UPP
- `aligned.fasta`: Final alignment
- `insertion_columns.txt`: Information about inserted columns
- `alignment_insertion_columns.txt`: Alignment with insertion annotations

### WITCH
- `aligned.fasta`: Final alignment
- `tree_decomp/`: Directory containing:
  * Cluster-specific alignments
  * HMM models
  * Decomposition tree
- `weights.txt`: HMM weights for each sequence

### CV-WITCH
- `aligned.fasta`: Final alignment
- `cv_results.txt`: Cross-validation performance metrics
- `alignment_weights.json`: CV-derived HMM weights
- `fold_*/`: Directory for each CV fold containing:
  * Training/validation splits
  * Fold-specific HMMs and alignments

## Analyzing Results

The benchmark plots show three key metrics:
1. SPFN (Sum-of-Pairs False Negative rate)
2. SPFP (Sum-of-Pairs False Positive rate)
3. FN (False Negative rate)

Lower values indicate better performance. The plots demonstrate:
- WITCH's general improvement over UPP in most cases
- CV-WITCH's slightly higher error rates but potential for better generalization
- Dataset-specific performance patterns