import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from protein_predictor import get_protein_embedding, compute_hydrophobicity
from protein_analysis import ProteinAnalyzer, prepare_visualization_data
import asyncio
import plotly.graph_objects as go
import plotly.express as px

st.set_page_config(
    page_title="Protein Analysis AI",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Protein Analysis AI")
st.markdown("""
### Analyze Your Protein Sequence
This tool helps you understand your protein's properties:
- 🌊 **Hydrophobicity**: How water-loving or water-fearing your protein is
- 🧪 **Function**: Predicts what your protein might do (e.g., membrane protein, enzyme)
- 🔍 **Domains**: Identifies functional regions in your protein
- 🧬 **Structure**: Estimates secondary structure elements

Try these example sequences:
```
Membrane protein: LLLLLLFFFFWWWWVVVVIII
DNA-binding:     KKKKRRRRHHHHDDDDEEEE
Insulin:         FVNQHLCGSHLVEAL
```
""")

# Initialize protein analyzer
if 'protein_analyzer' not in st.session_state:
    st.session_state.protein_analyzer = ProteinAnalyzer()

# Store sequences & analysis results in session state
if "sequence" not in st.session_state:
    st.session_state.sequence = ""
if "mutated_sequence" not in st.session_state:
    st.session_state.mutated_sequence = ""
if "analysis_results" not in st.session_state:
    st.session_state.analysis_results = None
if "mutated_analysis_results" not in st.session_state:
    st.session_state.mutated_analysis_results = None
if "has_prediction" not in st.session_state:
    st.session_state.has_prediction = False

st.session_state.sequence = st.text_area("Enter a protein sequence:", st.session_state.sequence)
if st.button("Predict"):
    sequence = st.session_state.sequence.upper().strip()
    if not sequence:
        st.error("⚠️ Please enter a protein sequence")
    else:
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        invalid_chars = [char for char in sequence if char not in valid_amino_acids]
        
        if invalid_chars:
            st.error(f"⚠️ Invalid characters detected: {', '.join(set(invalid_chars))}")
        else:
            with st.spinner("Analyzing protein sequence... This may take a few moments."):
                # Get hydrophobicity score
                hydrophobicity = compute_hydrophobicity(sequence)
                
                # Run real protein analysis
                analysis_results = asyncio.run(st.session_state.protein_analyzer.analyze_protein(sequence))
                st.session_state.analysis_results = analysis_results
                
                # Prepare data for visualization
                viz_data = prepare_visualization_data(analysis_results)
                
                # Store original analysis in session state for mutation comparison
                st.session_state.original_viz_data = viz_data
                
                # Display results
                st.subheader("Results:")
                
                # 1. Hydrophobicity Analysis
                st.write("### 🌊 Hydrophobicity Score")
                # Convert hydrophobicity from [-1, 1] to [0, 100]
                hydro_percent = min(max((hydrophobicity + 1) * 50, 0), 100)  # Clamp between 0 and 100
                st.write(f"Score: {hydro_percent:.1f}/100")
                st.write("- 0-49: Hydrophilic (water-loving)")
                st.write("- 50: Neutral")
                st.write("- 51-100: Hydrophobic (water-fearing)")
                
                # Create gradient from 0 to 100
                gradient = np.linspace(0, 100, 100)
                fig = go.Figure()
                
                # Add gradient bar
                fig.add_trace(go.Heatmap(
                    z=[gradient],
                    colorscale=[
                        [0, 'blue'],    # Hydrophilic
                        [0.5, 'white'], # Neutral
                        [1, 'red']      # Hydrophobic
                    ],
                    showscale=False
                ))
                
                # Add marker for the current sequence
                fig.add_trace(go.Scatter(
                    x=[hydro_percent],
                    y=[0],
                    mode='markers',
                    marker=dict(size=15, color='black'),
                    showlegend=False
                ))
                
                # Update layout
                fig.update_layout(
                    height=100,
                    margin=dict(l=50, r=50, t=0, b=0),
                    yaxis=dict(showticklabels=False),
                    xaxis=dict(
                        ticktext=['Hydrophilic', 'Neutral', 'Hydrophobic'],
                        tickvals=[25, 50, 75],
                        tickmode='array'
                    )
                )
                st.plotly_chart(fig)

                # 2. Protein Function Prediction
                st.write("### 🧪 Predicted Functions")
                if viz_data['function_scores']:
                    # Create bar chart
                    functions = list(viz_data['function_scores'].keys())
                    scores = list(viz_data['function_scores'].values())
                    
                    fig = go.Figure()
                    fig.add_trace(go.Bar(
                        x=functions,
                        y=scores,
                        marker_color='royalblue'
                    ))
                    
                    fig.update_layout(
                        height=300,
                        title='Function Predictions',
                        yaxis_title='Confidence (%)'
                    )
                    st.plotly_chart(fig)
                    
                    # Explain each predicted function
                    st.write("**What is the function of this protein?**")
                    for func, score in viz_data['function_scores'].items():
                        st.write(f"- This protein is likely a **{func.lower()}** ({score:.1f}% confidence)")
                        if func == "Membrane protein":
                            st.write("  - It contains many water-fearing amino acids, suggesting it spans or associates with cell membranes")
                        elif func == "DNA/RNA binding":
                            st.write("  - It has many charged amino acids, allowing it to interact with DNA or RNA molecules")
                        elif func == "Protein-protein interaction":
                            st.write("  - It has surface-exposed amino acids that can form bonds with other proteins")
                else:
                    st.write("No function predictions available")
                
                # 3. Domain Architecture
                st.write("### 🔍 Functional Domains")
                if viz_data['domain_architecture']:
                    for domain in viz_data['domain_architecture']:
                        st.write(f"- {domain}")
                else:
                    st.write("No clear domains detected in this sequence")

                # 4. Secondary Structure
                st.write("### 🧬 Secondary Structure")
                ss_percentages = viz_data['secondary_structure_percentages']
                if ss_percentages:
                    # Create pie chart
                    fig = px.pie(
                        values=list(ss_percentages.values()),
                        names=list(ss_percentages.keys()),
                        title='Predicted Secondary Structure',
                        color_discrete_map={
                            'Helix': '#FF9F1C',  # Orange
                            'Sheet': '#2EC4B6',  # Teal
                            'Coil': '#CCCCCC'    # Gray
                        }
                    )
                    fig.update_layout(height=400)
                    st.plotly_chart(fig)
                    
                    # Add explanation
                    st.write("**Structure Breakdown:**")
                    st.write(f"- Alpha Helix: {ss_percentages['Helix']:.1f}%")
                    st.write(f"- Beta Sheet: {ss_percentages['Sheet']:.1f}%")
                    st.write(f"- Random Coil: {ss_percentages['Coil']:.1f}%")
                else:
                    st.write("No secondary structure predictions available")

                # Add explanation of predictions
                st.write("\n### What do these predictions mean?")
                st.write("""
                - **Hydrophobicity**: Shows how water-loving (0-49) or water-fearing (51-100) the protein is. This affects where the protein is likely to be found in the cell.
                - **Functions**: Predicts potential roles based on amino acid patterns. Higher scores mean more confidence in the prediction.
                - **Domains**: Special regions in the protein that have specific roles, like binding to other molecules or spanning membranes.
                - **Structure**: Predicts how the protein folds (Helix = spiral, Sheet = flat, Coil = flexible).
                """)
                
                st.session_state.has_prediction = True

if st.session_state.has_prediction:
    st.subheader("Compare Mutations")
    st.write("Enter a mutated version of the sequence to compare properties:")
    st.session_state.mutated_sequence = st.text_input("Mutated sequence:", st.session_state.mutated_sequence).upper()
    
    if st.button("Compare Mutation"):
        mutated_sequence = st.session_state.mutated_sequence
        if not mutated_sequence:
            st.warning("⚠️ Please enter a mutated sequence to compare")
        elif len(mutated_sequence) != len(st.session_state.sequence):
            st.warning("⚠️ The mutated sequence must be the same length as the original")
        else:
            with st.spinner("Analyzing mutation effects..."):
                # Get original data from session state
                orig_viz_data = st.session_state.original_viz_data
                
                # Run analysis on mutated sequence
                mutated_results = asyncio.run(st.session_state.protein_analyzer.analyze_protein(mutated_sequence))
                mutated_viz_data = prepare_visualization_data(mutated_results)
                
                # Show comparisons
                st.subheader("🔄 Mutation Effects")
                
                # Compare hydrophobicity
                orig_hydro = compute_hydrophobicity(st.session_state.sequence)
                mut_hydro = compute_hydrophobicity(mutated_sequence)
                
                st.write("### 🌊 Hydrophobicity Change")
                orig_hydro_percent = min(max((orig_hydro + 1) * 50, 0), 100)  # Clamp between 0 and 100
                mut_hydro_percent = min(max((mut_hydro + 1) * 50, 0), 100)   # Clamp between 0 and 100
                hydro_change = mut_hydro_percent - orig_hydro_percent
                
                st.write(f"Original: {orig_hydro_percent:.1f}/100")
                st.write(f"Mutated: {mut_hydro_percent:.1f}/100")
                st.write(f"Change: {hydro_change:+.1f}")
                
                # Create comparison gradient
                fig = go.Figure()
                gradient = np.linspace(0, 100, 100)
                
                # Add gradient bar
                fig.add_trace(go.Heatmap(
                    z=[gradient],
                    colorscale=[
                        [0, 'blue'],    # Hydrophilic
                        [0.5, 'white'], # Neutral
                        [1, 'red']      # Hydrophobic
                    ],
                    showscale=False
                ))
                
                # Add markers for both sequences
                fig.add_trace(go.Scatter(
                    x=[orig_hydro_percent],
                    y=[0],
                    mode='markers',
                    name='Original',
                    marker=dict(size=15, color='black')
                ))
                
                fig.add_trace(go.Scatter(
                    x=[mut_hydro_percent],
                    y=[0],
                    mode='markers',
                    name='Mutated',
                    marker=dict(size=15, color='darkgreen')
                ))
                
                fig.update_layout(
                    height=100,
                    margin=dict(l=50, r=50, t=0, b=0),
                    yaxis=dict(showticklabels=False),
                    xaxis=dict(
                        ticktext=['Hydrophilic', 'Neutral', 'Hydrophobic'],
                        tickvals=[25, 50, 75],
                        tickmode='array'
                    ),
                    showlegend=True,
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.02,
                        xanchor="right",
                        x=1
                    )
                )
                st.plotly_chart(fig)
                
                if abs(hydro_change) > 25:
                    st.warning("⚠️ Large change in hydrophobicity may affect protein stability")
                
                # Compare functions
                st.write("### 🧪 Function Changes")
                orig_funcs = orig_viz_data['function_scores']
                mut_funcs = mutated_viz_data['function_scores']
                
                if orig_funcs or mut_funcs:  # Show if either has predictions
                    # Create comparison bar chart
                    functions = list(set(orig_funcs.keys()) | set(mut_funcs.keys()))
                    orig_scores = [orig_funcs.get(f, 0) for f in functions]
                    mut_scores = [mut_funcs.get(f, 0) for f in functions]
                    
                    fig = go.Figure()
                    fig.add_trace(go.Bar(
                        name='Original',
                        x=functions,
                        y=orig_scores,
                        marker_color='royalblue'
                    ))
                    fig.add_trace(go.Bar(
                        name='Mutated',
                        x=functions,
                        y=mut_scores,
                        marker_color='darkgreen'
                    ))
                    
                    fig.update_layout(
                        barmode='group',
                        height=300,
                        title='Function Prediction Comparison',
                        yaxis_title='Confidence (%)'
                    )
                    st.plotly_chart(fig)
                    
                    # Show significant changes
                    significant_changes = []
                    for func in functions:
                        orig = orig_funcs.get(func, 0)
                        mut = mut_funcs.get(func, 0)
                        if abs(mut - orig) > 10:
                            significant_changes.append(f"- {func}: {orig:.1f}% → {mut:.1f}% ({mut-orig:+.1f}%)")
                    
                    if significant_changes:
                        st.write("**Significant Changes:**")
                        for change in significant_changes:
                            st.write(change)
                    else:
                        st.write("No significant changes in function predictions")
                else:
                    st.write("No function predictions available")
                
                # Compare domains
                st.write("### 🔍 Domain Changes")
                orig_domains = orig_viz_data.get('domain_architecture', [])
                mut_domains = mutated_viz_data.get('domain_architecture', [])
                
                col1, col2 = st.columns(2)
                with col1:
                    st.write("**Original Domains:**")
                    if orig_domains:
                        for domain in orig_domains:
                            st.write(f"- {domain}")
                    else:
                        st.write("No domains detected")
                
                with col2:
                    st.write("**Mutated Domains:**")
                    if mut_domains:
                        for domain in mut_domains:
                            st.write(f"- {domain}")
                    else:
                        st.write("No domains detected")
                
                # Compare secondary structure
                st.write("### 🧬 Structure Changes")
                orig_ss = orig_viz_data.get('secondary_structure_percentages', {})
                mut_ss = mutated_viz_data.get('secondary_structure_percentages', {})
                
                if orig_ss and mut_ss:
                    # Create side-by-side pie charts
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        fig1 = px.pie(
                            values=list(orig_ss.values()),
                            names=list(orig_ss.keys()),
                            title='Original Structure',
                            color_discrete_map={
                                'Helix': '#FF9F1C',
                                'Sheet': '#2EC4B6',
                                'Coil': '#CCCCCC'
                            }
                        )
                        fig1.update_layout(height=300)
                        st.plotly_chart(fig1)
                    
                    with col2:
                        fig2 = px.pie(
                            values=list(mut_ss.values()),
                            names=list(mut_ss.keys()),
                            title='Mutated Structure',
                            color_discrete_map={
                                'Helix': '#FF9F1C',
                                'Sheet': '#2EC4B6',
                                'Coil': '#CCCCCC'
                            }
                        )
                        fig2.update_layout(height=300)
                        st.plotly_chart(fig2)
                    
                    # Show significant changes
                    changes = []
                    for struct in ['Helix', 'Sheet', 'Coil']:
                        if struct in orig_ss and struct in mut_ss:
                            change = mut_ss[struct] - orig_ss[struct]
                            if abs(change) > 5:  # Only show significant changes
                                changes.append(f"- {struct}: {orig_ss[struct]:.1f}% → {mut_ss[struct]:.1f}% ({change:+.1f}%)")
                    
                    if changes:
                        st.write("**Significant Changes:**")
                        for change in changes:
                            st.write(change)
                    else:
                        st.write("No significant changes in structure")
                else:
                    st.write("No structure predictions available")
