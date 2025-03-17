import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from protein_predictor import get_protein_embedding, compute_hydrophobicity

def softmax(x):
    """Compute softmax values for an array."""
    exp_x = np.exp(x - np.max(x))  # Subtract max for numerical stability
    return exp_x / exp_x.sum()

st.title("🧬 Protein Property Prediction AI")

# Store sequences & hydrophobicity in session state
if "sequence" not in st.session_state:
    st.session_state.sequence = ""
if "mutated_sequence" not in st.session_state:
    st.session_state.mutated_sequence = ""
if "hydrophobicity" not in st.session_state:
    st.session_state.hydrophobicity = None
if "mutated_hydrophobicity" not in st.session_state:
    st.session_state.mutated_hydrophobicity = None

# **Predefined Protein Function Embeddings**
protein_function_db = {
    "enzyme": np.array([0.8, 0.2, 0.5, 0.6, 0.9] * 10),
    "transport": np.array([0.5, 0.7, 0.4, 0.3, 0.6] * 10),
    "signaling": np.array([0.2, 0.4, 0.8, 0.9, 0.3] * 10),
}

def predict_protein_function(embedding):
    embedding = embedding[:50]  
    similarities = {func: np.dot(embedding, vec) for func, vec in protein_function_db.items()}
    confidence_scores = softmax(np.array(list(similarities.values()))) * 100  
    predicted_function = max(similarities, key=similarities.get)
    predicted_confidence = confidence_scores[list(similarities.keys()).index(predicted_function)]
    return predicted_function, predicted_confidence

def predict_secondary_structure(sequence):
    helix_count = sum(sequence.count(aa) for aa in "ALMQE")  
    sheet_count = sum(sequence.count(aa) for aa in "VIFWY")  
    total_count = len(sequence)
    helix_fraction = helix_count / total_count
    sheet_fraction = sheet_count / total_count
    coil_fraction = 1 - (helix_fraction + sheet_fraction)
    return helix_fraction, sheet_fraction, coil_fraction

def predict_functional_domain(sequence):
    transmembrane_score = sum(sequence.count(aa) for aa in "LVIFW") / len(sequence)
    catalytic_score = sum(sequence.count(aa) for aa in "DEHRK") / len(sequence)
    signal_score = sum(sequence.count(aa) for aa in "ACG") / len(sequence)
    
    scores = {
        "transmembrane domain": transmembrane_score,
        "catalytic domain": catalytic_score,
        "signal peptide": signal_score,
    }
    
    best_domain = max(scores, key=scores.get)
    best_score = scores[best_domain]
    
    if best_score < 0.1:
        return "no clear domain", "meaning its function is uncertain."
    else:
        return best_domain, f"suggesting it could be involved in {best_domain} activity."

st.session_state.sequence = st.text_area("Enter a protein sequence:", st.session_state.sequence)
st.session_state.mutated_sequence = st.text_input("Enter a mutated version of your sequence:", st.session_state.mutated_sequence)

if st.button("Predict"):
    sequence = st.session_state.sequence.upper()
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_chars = [char for char in sequence if char not in valid_amino_acids]
    
    if invalid_chars:
        st.error(f"⚠️ Invalid characters detected: {', '.join(set(invalid_chars))}")
    else:
        st.subheader("Results:")
        st.session_state.hydrophobicity = compute_hydrophobicity(sequence)
        st.write(f"🔹 **Hydrophobicity Score:** {st.session_state.hydrophobicity:.2f}")

        # **Hydrophobicity Scale**
        gradient = np.linspace(-4, 4, 300)
        cmap = plt.cm.coolwarm  
        st.subheader("🌊 Hydrophobicity Scale")
        fig, ax = plt.subplots(figsize=(8, 1))
        ax.imshow(gradient.reshape(1, -1), aspect="auto", cmap=cmap, extent=[-4, 4, -0.3, 0.3])
        ax.scatter(st.session_state.hydrophobicity, 0, color="black", edgecolors="white", s=150, zorder=3)
        ax.set_xlim(-4, 4)
        ax.axvline(0, color="black", linewidth=1, linestyle="dashed")
        ax.set_yticks([])  # Remove Y-Axis Labels
        st.pyplot(fig)

        # **AI-Predicted Function**
        st.subheader("🧪 AI-Predicted Protein Function")
        embedding = get_protein_embedding(sequence).numpy().flatten()
        predicted_function, confidence = predict_protein_function(embedding)
        st.write(f"🔹 **This protein is {confidence:.1f}% similar to known {predicted_function} proteins.**")

        # **Functional Domain Prediction**
        st.subheader("🔍 Functional Domain Prediction")
        domain, explanation = predict_functional_domain(sequence)
        st.write(f"🔹 **Based on its structure, this protein may contain a {domain}, {explanation}**")

        # **Secondary Structure**
        st.subheader("🧬 Predicted Secondary Structure")
        helix_fraction, sheet_fraction, coil_fraction = predict_secondary_structure(sequence)
        fig, ax = plt.subplots()
        ax.pie([helix_fraction, sheet_fraction, coil_fraction], labels=["Helix", "Sheet", "Coil"],
               autopct="%1.1f%%", colors=["orange", "blue", "gray"])
        st.pyplot(fig)

if st.button("Mutate"):
    mutated_sequence = st.session_state.mutated_sequence.upper()
    if len(mutated_sequence) != len(st.session_state.sequence):
        st.warning("⚠️ The mutated sequence must be the same length as the original.")
    else:
        st.subheader("🔄 Mutation Results:")
        
        st.session_state.mutated_hydrophobicity = compute_hydrophobicity(mutated_sequence)
        st.write(f"🔹 **Mutated Hydrophobicity Score:** {st.session_state.mutated_hydrophobicity:.2f}")

        # **Mutated Hydrophobicity Scale**
        st.subheader("🌊 Mutated Hydrophobicity Scale")
        fig, ax = plt.subplots(figsize=(8, 1))
        ax.imshow(gradient.reshape(1, -1), aspect="auto", cmap=cmap, extent=[-4, 4, -0.3, 0.3])
        ax.scatter(st.session_state.mutated_hydrophobicity, 0, color="red", edgecolors="white", s=150, zorder=3)
        ax.set_xlim(-4, 4)
        ax.axvline(0, color="black", linewidth=1, linestyle="dashed")
        ax.set_yticks([])
        st.pyplot(fig)

        # **AI-Predicted Function**
        st.subheader("🧪 Mutated AI-Predicted Protein Function")
        mutated_embedding = get_protein_embedding(mutated_sequence).numpy().flatten()
        mutated_function, mutated_confidence = predict_protein_function(mutated_embedding)
        st.write(f"🔹 **This mutated protein is {mutated_confidence:.1f}% similar to known {mutated_function} proteins.**")

        # **Mutated Secondary Structure**
        st.subheader("🧬 Mutated Secondary Structure")
        mutated_helix, mutated_sheet, mutated_coil = predict_secondary_structure(mutated_sequence)
        fig, ax = plt.subplots()
        ax.pie([mutated_helix, mutated_sheet, mutated_coil], labels=["Helix", "Sheet", "Coil"],
               autopct="%1.1f%%", colors=["orange", "blue", "gray"])
        st.pyplot(fig)


