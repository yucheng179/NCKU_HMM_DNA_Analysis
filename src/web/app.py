from flask import Flask, render_template, request, jsonify
import subprocess

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/run_model', methods=['POST'])
def run_model():
    model = request.json['model']
    cmd_map = {
    'zero': './models/zero_order.exe',
    'first': './models/first_order.exe',
    'second': './models/second_order.exe',
    'hmm': './models/hmm_project.exe',
    }

    if model not in cmd_map:
        return jsonify({'error': 'Invalid model name'}), 400

    try:
        result = subprocess.check_output(cmd_map[model], stderr=subprocess.STDOUT, timeout=60)
        return jsonify({'output': result.decode('utf-8')})
    except subprocess.CalledProcessError as e:
        return jsonify({'error': e.output.decode('utf-8')}), 500

if __name__ == '__main__':
    app.run(debug=True)
