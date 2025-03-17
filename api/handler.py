from flask import Flask, request, jsonify
import os
import sys

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from protein_analysis import ProteinAnalyzer

app = Flask(__name__)
analyzer = ProteinAnalyzer()

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
