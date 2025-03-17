from flask import Flask, request, jsonify

app = Flask(__name__)

@app.route('/', methods=['GET'])
def home():
    return jsonify({"message": "Protein Analysis API is running"})

@app.route('/api/test', methods=['GET'])
def test():
    return jsonify({"status": "ok"})

@app.route('/api/analyze', methods=['POST'])
def analyze():
    try:
        data = request.get_json()
        if not data or 'sequence' not in data:
            return jsonify({'error': 'No sequence provided'}), 400
        
        sequence = data['sequence'].upper()
        return jsonify({
            'sequence': sequence,
            'length': len(sequence)
        })
    except Exception as e:
        print(f"Error: {str(e)}")
        return jsonify({'error': str(e)}), 500
