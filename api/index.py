from flask import Flask, request, jsonify
import sys
import os

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from protein_analysis import ProteinAnalyzer

app = Flask(__name__)
analyzer = ProteinAnalyzer()

@app.route('/', methods=['GET'])
def home():
    return """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Protein Analysis API</title>
        <style>
            body { font-family: Arial, sans-serif; max-width: 800px; margin: 0 auto; padding: 20px; }
            pre { background: #f5f5f5; padding: 15px; border-radius: 5px; }
        </style>
    </head>
    <body>
        <h1>🧬 Protein Analysis API</h1>
        <p>Use the /api/analyze endpoint to analyze protein sequences.</p>
        <h2>Example Usage:</h2>
        <pre>
POST /api/analyze
Content-Type: application/json

{
    "sequence": "FVNQHLCGSHLVEAL"
}
        </pre>
    </body>
    </html>
    """

@app.route('/api/analyze', methods=['POST'])
def analyze():
    data = request.get_json()
    if not data or 'sequence' not in data:
        return jsonify({'error': 'No sequence provided'}), 400
    
    sequence = data['sequence'].upper()
    
    # Validate sequence
    if not all(aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in sequence):
        return jsonify({'error': 'Invalid sequence. Use only valid amino acid letters.'}), 400
    
    try:
        # Analyze sequence
        results = analyzer.analyze_protein(sequence)
        return jsonify(results)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run()
