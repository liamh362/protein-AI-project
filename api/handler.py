from flask import Flask, request, jsonify
import numpy as np

app = Flask(__name__)

class SimpleProteinAnalyzer:
    def __init__(self):
        self.hydrophobic = set("VILMFYW")
        self.polar = set("STNQ")
        self.charged = set("DEKR")
        
    def analyze_protein(self, sequence: str) -> dict:
        """Basic protein sequence analysis."""
        total_len = len(sequence)
        
        # Calculate amino acid composition
        hydrophobic_content = sum(aa in self.hydrophobic for aa in sequence) / total_len
        polar_content = sum(aa in self.polar for aa in sequence) / total_len
        charged_content = sum(aa in self.charged for aa in sequence) / total_len
        
        # Calculate molecular weight (approximate)
        avg_aa_weight = 110  # Average amino acid weight in Daltons
        molecular_weight = total_len * avg_aa_weight
        
        return {
            'sequence_length': total_len,
            'molecular_weight': molecular_weight,
            'composition': {
                'hydrophobic': round(hydrophobic_content * 100, 2),
                'polar': round(polar_content * 100, 2),
                'charged': round(charged_content * 100, 2)
            },
            'sequence': sequence
        }

analyzer = SimpleProteinAnalyzer()

@app.route('/')
def home():
    return jsonify({"message": "Protein Analysis API is running"})

@app.route('/api/analyze', methods=['POST'])
def analyze():
    try:
        data = request.get_json()
        if not data or 'sequence' not in data:
            return jsonify({'error': 'No sequence provided'}), 400
        
        sequence = data['sequence'].upper()
        
        # Validate sequence
        if not all(aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in sequence):
            return jsonify({'error': 'Invalid sequence. Use only valid amino acid letters.'}), 400
        
        # Analyze sequence
        results = analyzer.analyze_protein(sequence)
        return jsonify(results)
    except Exception as e:
        print(f"Error: {str(e)}")
        return jsonify({'error': str(e)}), 500

app.debug = True
