# RNA Half-Life Prediction

Predicting mRNA half-life and determining the effect of synonymous mutations on half-life predictions using deep learning.

## Overview

This project uses a dual-input convolutional neural network (CNN) to predict mRNA half-life based on sequence features and RNA secondary structure. The model combines:
- **Sequence features**: One-hot encoded RNA sequences with base-pairing probabilities
- **Biological features**: 13 computed features including GC content, AU-rich elements (ARE), codon adaptation index (CAI), miRNA binding sites, and ORF characteristics

The project also explores how synonymous mutations affect predicted mRNA stability using *CASP9* (Caspase 9) as a test case.

## Requirements

### R Packages
```r
# Data manipulation
library(tidyverse)
library(readr)

# Bioinformatics
library(Biostrings)
library(seqinr)

# Machine learning
library(keras3)
library(abind)
library(rsample)
```

### External Tools
- **ViennaRNA** (RNAfold): Required for calculating RNA secondary structure and base-pairing probabilities
  - Install from: https://www.tbi.univie.ac.at/RNA/

### Data Files (not included)
- `rna_half_lives.csv`: mRNA half-life measurements across multiple samples
- `all_human_rna.fna`: FASTA file containing human RNA sequences

## Pipeline Overview

The analysis consists of 8 sequential scripts:

### 1. Dataset Preprocessing (`01_dataset_preprocessing.R`)
- Loads mRNA half-life data and sequences
- Filters genes with available sequences
- Converts sequences to RNA format (T → U)
- Calculates mean half-life across samples
- Generates input files for ViennaRNA structure prediction
- Combines ViennaRNA output with sequence data

**Output**: `cnn_input_with_pairing.csv`, `cnn_input_grouped.rds`

### 2. Feature Engineering (`02_prepare_data.R`)
- Filters sequences by length (500-5000 bp)
- Calculates 13 biological features:
  - Sequence composition (GC%, AU%)
  - 3' UTR features (GC%, pairing probability)
  - AU-rich element (ARE) score
  - Poly(A) signal presence
  - ORF detection and characteristics
  - Codon Adaptation Index (CAI)
  - miRNA binding site density (top 20 human miRNAs)
- One-hot encodes sequences
- Pads sequences to uniform length (3500 bp)
- Splits data into train/validation/test sets (70/15/15)
- Normalizes features and targets (log transformation)

**Output**: `cnn_training_data_with_features.rds`

### 3. Model Training (`03_cnn_training.R`)
- Builds dual-input CNN architecture:
  - **Sequence pathway**: Dilated 1D convolutions (rates: 1, 2, 4) with global max pooling
  - **Feature pathway**: Dense layers with batch normalization
  - Combined pathways for final prediction
- Uses dropout for regularization
- Implements early stopping and learning rate reduction
- Evaluates on test set with correlation metrics

**Output**: `rna_halflife_improved_weights.weights.h5`, `improved_model_predictions.rds`

### 4. Generate CASP9 Mutants (`04_mutate_casp9.R`)
- Selects *CASP9* (Caspase 9) as test gene
- Identifies coding sequence (CDS) via ORF detection
- Creates 75 synonymous mutants in 300 bp region:
  - 25 mutants with 10 mutations each
  - 25 mutants with 50 mutations each
  - 25 mutants with 100 mutations each
- Maintains amino acid sequence (silent mutations only)
- Generates FASTA file for ViennaRNA analysis

**Output**: `caspase9_mutants.rds`, `caspase9_all_sequences.fa`

### 5. Preprocess CASP9 Mutants (`05_casp9_preprocessing.R`)
- Loads ViennaRNA output for all mutants
- Calculates all 13 features for wildtype + mutants
- Normalizes features using training data parameters
- Prepares sequence tensors with padding and masking

**Output**: `caspase9_prepared_for_prediction.rds`

### 6. Predict CASP9 Half-Lives (`06_casp9_mut_testing.R`)
- Loads trained model weights
- Makes predictions on wildtype and all mutants
- Transforms predictions back to original scale
- Calculates change from wildtype for each mutant
- Performs diagnostic checks on prediction consistency

**Output**: `caspase9_predictions.rds`, `caspase9_predictions.csv`

### 7. Statistical Analysis (`07_analysis.R`)
- Summarizes mutation effects by mutation count
- Identifies features that correlate with stability changes
- Performs regression analysis (linear and quadratic models)
- Conducts statistical tests (t-tests, ANOVA)
- Compares mutation-sensitive features vs. overall feature importance
- Identifies extreme mutants (most/least stable)

**Output**: `mutation_effects_for_plotting.csv`, `feature_correlations.rds`

### 8. Generate Figures (`08_make_plots.R`)
- **Figure 1**: Model performance (actual vs. predicted half-life)
- **Figure 2**: Feature importance in test set
- **Figure 3**: CASP9 mutation effects by mutation count
- **Figure 4**: Feature correlations with CASP9 stability changes

**Output**: PNG files in `Plots/` directory

## Usage

### Full Pipeline
```bash
# 1. Preprocess data
Rscript Scripts/01_dataset_preprocessing.R

# 2. Run ViennaRNA on filtered genes
RNAfold -p --noPS < filtered_genes.fa
# (Process output files into pairing_probs/ directory)

# 3. Prepare training data
Rscript Scripts/02_prepare_data.R

# 4. Train model
Rscript Scripts/03_cnn_training.R

# 5. Generate CASP9 mutants
Rscript Scripts/04_mutate_casp9.R

# 6. Run ViennaRNA on mutants
RNAfold -p --noPS < caspase9_all_sequences.fa
# (Process output files into pairing_probs_mutants/ directory)

# 7. Preprocess mutants
Rscript Scripts/05_casp9_preprocessing.R

# 8. Predict mutant half-lives
Rscript Scripts/06_casp9_mut_testing.R

# 9. Analyze results
Rscript Scripts/07_analysis.R

# 10. Generate plots
Rscript Scripts/08_make_plots.R
```

## Model Architecture

### Sequence Input Branch
- Input: (3500, 6) - [A, C, G, U, pairing_prob, mask]
- Conv1D (64 filters, kernel=3, dilation=1)
- Conv1D (64 filters, kernel=3, dilation=2)
- Conv1D (128 filters, kernel=3, dilation=4)
- Global Max Pooling

### Feature Input Branch
- Input: (13,) - biological features
- Dense (128) + BatchNorm + Dropout (0.3)
- Dense (64) + BatchNorm + Dropout (0.3)
- Dense (32) + Dropout (0.2)

### Combined Output
- Concatenate both branches
- Dense (128) + Dropout (0.5)
- Dense (64) + Dropout (0.4)
- Dense (1) - final half-life prediction

## Key Features

### Biological Features (n=13)
1. **Sequence length**: Total nucleotides
2. **GC content**: Overall G+C percentage
3. **AU content**: Overall A+U percentage
4. **ARE score**: AUUUA motif density in 3' UTR
5. **Poly(A) signal**: Presence of canonical polyadenylation signals
6. **3' UTR GC content**: GC% in last 500 bp
7. **Mean pairing probability**: Average base-pairing probability
8. **3' UTR pairing**: Mean pairing probability in 3' UTR
9. **ORF length**: Longest open reading frame
10. **ORF ratio**: ORF length / total length
11. **CAI**: Codon Adaptation Index (translation efficiency)
12. **miRNA sites**: Number of miRNA seed matches (top 20 human miRNAs)
13. **miRNA density**: miRNA sites per 1000 bp

## Results

### Model Performance (Test Set)
- Pearson correlation: r = 0.227
- Spearman correlation: ρ = [value from your results]
- Median absolute error: [value] hours

### CASP9 Mutation Analysis
- Average effect per synonymous mutation: ~0.0001-0.001 hours
- Maximum observed change: [value]% from wildtype
- Correlation between mutation count and stability: r = [value]

### Key Findings
- Long-lived mRNAs (>15 hrs) show [specific behavior]
- Top features driving predictions: [list from Figure 2]
- Synonymous mutations have minimal effect on predicted stability
- Features most sensitive to mutations: CAI, miRNA binding sites

## Directory Structure
```
rna-halflife-prediction/
├── Scripts/
│   ├── 01_dataset_preprocessing.R
│   ├── 02_prepare_data.R
│   ├── 03_cnn_training.R
│   ├── 04_mutate_casp9.R
│   ├── 05_casp9_preprocessing.R
│   ├── 06_casp9_mut_testing.R
│   ├── 07_analysis.R
│   └── 08_make_plots.R
├── Plots/
│   ├── figure1.png
│   ├── figure2.png
│   ├── figure3.png
│   └── figure4.png
├── Writing/
└── README.md
```

## Contact

Jason Summers - [jasonsummers012]

## Acknowledgments

- ViennaRNA package for RNA secondary structure prediction
- Keras/TensorFlow for deep learning framework
- Biostrings and seqinr for sequence analysis
