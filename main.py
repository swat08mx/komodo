from flask import Flask, render_template, request, redirect, url_for
import os
import subprocess
from werkzeug.utils import secure_filename


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads/'

PIPELINES = {
    'RNA-Seq': 'pipelines/single_cell.py',
    'Variant Calling': 'pipelines/variant_calling.py',
    'Metagenomics': 'pipelines/metagenomics.py'
}

@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    if request.method == 'POST':
        pipeline_name = request.form.get('pipeline')
        file = request.files.get('sample_file')
        if file and pipeline_name:
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)
            script_path = PIPELINES.get(pipeline_name)
            if script_path:
                try:
                    # Run the Python script and capture output
                    completed = subprocess.run(
                        ['python', script_path, filepath],
                        capture_output=True, text=True, check=True
                    )
                    result = completed.stdout
                except subprocess.CalledProcessError as e:
                    result = f"Error: {e.stderr}"
            else:
                result = "Invalid pipeline selected."
    return render_template('index.html', pipelines=PIPELINES.keys(), result=result)

@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        file = request.files['file']
        if file:
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return redirect(url_for('download_file', name=filename))
    return render_template('index.html')


if __name__ == '__main__':
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    app.run(host='0.0.0.0', port='5000', debug=True)