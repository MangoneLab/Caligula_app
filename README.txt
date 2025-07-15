An AI-Guided Framework Reveals Conserved Rules Governing microRNA Strand Selection
Dalton Meadows1,2, 
Hailee Hargis1,2,
Amanda Ellis1,2,
Heewook Lee1,
Marco Mangone1,3

Author affiliations
1) The Biodesign Institute at Arizona State University, 1001 S McAllister Ave, Tempe, AZ
2) School of Life Sciences, Arizona State University, 427 E Tyler Mall, Tempe, AZ
3) Corresponding Author

Requirements:
Python 3.10
pandas	1.5.3	Reading/writing CSVs (batch.py)
joblib	1.2.0	Loading pre-trained model (model.pkl)
numpy	1.23.5	Core numerical operations
scikit-learn	1.2.2	Model training and prediction


Two python scripts are provided: caligula_single.py and caliguyla_batch.py

##################
caligula_single.py
##################

usage: python caligula_single.py

This script is a standalone terminal application that predicts which strand (5p or 3p) of a given miRNA hairpin is more likely to be selected for RISC loading, based solely on structural and sequence-derived features from the uppercase-labeled regions.

Core Workflow:

User Input:
Prompts the user to enter a hairpin sequence in which the two uppercase regions represent candidate 5p and 3p strands. (For simplicity we hard-coded let-7 miRNA sequence and its predicted hairpin structure).

Feature Extraction:
Internally calls extract_strand_features.py to compute ~77 biologically informed features derived only from the uppercase regions. These include:
GC content, nucleotide frequencies
Positional biases (e.g., seed pairing, 5′ U)
Local structure features (e.g., arm pairing probability, loop distance)
Thermodynamic features (optionally simplified versions without ViennaRNA dependency)

Model Prediction:
Loads a pre-trained machine learning model (model.pkl)
Applies the model to the extracted features
Outputs a prediction: 5p or 3p, along with a confidence score.


##################
caligula_batch.py
##################

Usage: python caligula_batch.py input.csv

This script performs strand usage prediction (5p vs 3p) for multiple miRNA hairpins provided in a CSV file. It systematically processes each input, extracts features, and applies the trained Caligula model to predict the dominant strand at scale.

Input Format:
CSV file (passed as input.csv)
Must contain at least three columns:
name: name of the miRNA
sequence: the miRNA hairpin with exactly two uppercase regions, representing candidate 5p and 3p strands.
hairpin_structure: the mold predicted hairpin structure

Core Workflow:

Read Input File:
Loads the input CSV using pandas.

Feature Extraction Loop:
- Iterates over each hairpin sequence.
- For each, calls the feature extractor (e.g. extract_strand_features from extract_strand_features_fixed.py).
- Skips sequences that are malformed or missing two uppercase regions.

Model Loading and Prediction:
Loads the pre-trained machine learning model (model.pkl) using joblib.
Applies the model to the batch of extracted features.

Output Generation:
Appends predicted strand (5p or 3p) and confidence score to the original data.
Writes the annotated predictions to predictions.csv.
