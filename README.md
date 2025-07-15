# Caligula miRNA Strand Predictor
MangoneLab - Arizona State University

This script is a standalone terminal application that predicts which strand (5p or 3p) of a given miRNA hairpin is more likely to be selected for RISC loading, based solely on structural and sequence-derived features from the uppercase-labeled regions.

## Getting Started

1. Click **Code > Open with Codespaces > New codespace**
2. Wait for setup (~1 min)
3. In the terminal, run:

```bash
python caligula_single.py
```

4. Follow the prompt to enter a hairpin sequence and structure.


Core Workflow:

User Input:
Prompts the user to enter a hairpin sequence in which the two uppercase regions represent candidate 5p and 3p strands. For example, use seq:ugggaUGAGGUAGUAGGUUGUAUAGUUuuuCUAUACAAUCUACUGUCUUUCcua and structure: (((((((.((((((((((((((((.....))))))))))))))).))))))).

Feature Extraction:
Internally the script calls extract_strand_features.py to compute ~77 biologically informed features derived only from the uppercase regions. These include:
GC content, nucleotide frequencies
Positional biases (e.g., seed pairing, 5′ U)
Local structure features (e.g., arm pairing probability, loop distance)
Thermodynamic features (optionally simplified versions without ViennaRNA dependency)

Model Prediction:
Loads a pre-trained machine learning model (model.pkl)
Applies the model to the extracted features
Outputs a prediction: 5p or 3p, along with a confidence score.


