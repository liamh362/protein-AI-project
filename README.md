# 🧬 Protein Property Prediction AI

An interactive web application that analyzes protein sequences and predicts various properties using AI and bioinformatics techniques.

## Features

- 🔍 Protein Function Prediction
- 🌊 Hydrophobicity Analysis
- 🧪 Secondary Structure Prediction
- 🔄 Mutation Analysis
- 🎯 Functional Domain Identification

## Getting Started

### Prerequisites

- Python 3.x
- Streamlit
- NumPy
- Matplotlib

### Installation

1. Clone the repository:
```bash
git clone https://github.com/liamh362/protein-AI-project.git
cd protein-AI-project
```

2. Install dependencies:
```bash
pip install streamlit numpy matplotlib
```

3. Run the application:
```bash
streamlit run app.py
```

## Usage

1. Enter a protein sequence using standard amino acid codes (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
2. Click "Predict" to analyze the sequence
3. Optionally, enter a mutated version of your sequence to compare properties

## Features in Detail

### Protein Function Prediction
- Uses AI embeddings to compare against known protein functions
- Provides confidence scores for predictions

### Hydrophobicity Analysis
- Calculates overall hydrophobicity score
- Visual representation on a hydrophobicity scale

### Secondary Structure Analysis
- Predicts proportions of:
  - Alpha helices
  - Beta sheets
  - Random coils
- Visualizes structure distribution in a pie chart

### Mutation Analysis
- Compare properties between original and mutated sequences
- Analyze changes in:
  - Hydrophobicity
  - Predicted function
  - Secondary structure

## Contributing

Feel free to open issues or submit pull requests to improve the project.
