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
    return """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Protein Analysis Web App</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                max-width: 800px;
                margin: 0 auto;
                padding: 20px;
                line-height: 1.6;
            }
            .container {
                background-color: #f5f5f5;
                padding: 20px;
                border-radius: 8px;
                margin-top: 20px;
            }
            input[type="text"] {
                width: 100%;
                padding: 10px;
                margin: 10px 0;
                border: 1px solid #ddd;
                border-radius: 4px;
            }
            button {
                background-color: #4CAF50;
                color: white;
                padding: 10px 20px;
                border: none;
                border-radius: 4px;
                cursor: pointer;
            }
            button:hover {
                background-color: #45a049;
            }
            #result {
                margin-top: 20px;
                white-space: pre-wrap;
            }
            .error {
                color: red;
                margin-top: 10px;
            }
        </style>
    </head>
    <body>
        <h1>🧬 Protein Sequence Analyzer</h1>
        
        <div class="container">
            <h2>Enter Protein Sequence</h2>
            <input type="text" id="sequence" placeholder="Enter protein sequence (e.g., FVNQHLCGSHLVEAL)">
            <button onclick="analyzeProtein()">Analyze</button>
            <div id="result"></div>
        </div>

        <script>
        function analyzeProtein() {
            const sequence = document.getElementById('sequence').value;
            if (!sequence) {
                document.getElementById('result').innerHTML = '<div class="error">Please enter a protein sequence</div>';
                return;
            }

            fetch('/api/analyze', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({ sequence: sequence })
            })
            .then(response => response.json())
            .then(data => {
                if (data.error) {
                    document.getElementById('result').innerHTML = `<div class="error">${data.error}</div>`;
                } else {
                    const result = `
                        <h3>Analysis Results:</h3>
                        <p>Sequence: ${data.sequence}</p>
                        <p>Length: ${data.sequence_length} amino acids</p>
                        <p>Molecular Weight: ${data.molecular_weight} Daltons</p>
                        <h4>Composition:</h4>
                        <ul>
                            <li>Hydrophobic: ${data.composition.hydrophobic}%</li>
                            <li>Polar: ${data.composition.polar}%</li>
                            <li>Charged: ${data.composition.charged}%</li>
                        </ul>
                    `;
                    document.getElementById('result').innerHTML = result;
                }
            })
            .catch(error => {
                document.getElementById('result').innerHTML = `<div class="error">Error: ${error.message}</div>`;
            });
        }
        </script>
    </body>
    </html>
    """

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
