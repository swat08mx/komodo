from flask import Flask, render_template, request, redirect, url_for
import os

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads/'

# Example list of pipelines
PIPELINES = ['RNA-Seq', 'Variant Calling', 'Metagenomics']

@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    if request.method == 'POST':
        pipeline = request.form.get('pipeline')
        file = request.files.get('sample_file')
        if file and pipeline:
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)
            # Here, insert your backend pipeline execution logic
            # For demonstration, we mock a result:
            result = f"Pipeline '{pipeline}' executed on file '{file.filename}'."
    return render_template('index.html', pipelines=PIPELINES, result=result)

if __name__ == '__main__':
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    app.run(host='0.0.0.0', port='5000', debug=True)