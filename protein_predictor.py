from transformers import EsmModel, EsmTokenizer
import torch
import numpy as np

# Load the ESM-2 model
model_name = "facebook/esm2_t6_8M_UR50D"
tokenizer = EsmTokenizer.from_pretrained(model_name)
model = EsmModel.from_pretrained(model_name)

def get_protein_embedding(sequence):
    """Convert a protein sequence into an AI-generated numerical representation."""
    inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
    with torch.no_grad():
        outputs = model(**inputs)
    return outputs.last_hidden_state.mean(dim=1)  # Returns the AI-based feature embedding

# Define a hydrophobicity scale for each amino acid
hydrophobicity_dict = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
    'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

def compute_hydrophobicity(sequence):
    """Compute the average hydrophobicity of a protein sequence."""
    scores = [hydrophobicity_dict.get(aa, 0) for aa in sequence]  # Default 0 if unknown amino acid
    return np.mean(scores) if scores else 0

