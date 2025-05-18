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
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif', 'csv', 'tsv'}

# Create upload folder if it doesn't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def index():
    return render_template('index.html', pipelines=PIPELINES)


@app.route('/submission', methods=['POST'])
def submission():
    # Check if the post request has the file part
    if 'file' not in request.files:
        return 'No file part in the request'

    file = request.files['file']

    # If user does not select file
    if file.filename == '':
        return 'No selected file'

    if file and allowed_file(file.filename):
        # Secure the filename before storing
        filename = secure_filename(file.filename)
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        return 'File successfully uploaded, thank you!'

    return "INvalid file type"


if __name__ == '__main__':
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    app.run(host='0.0.0.0', port='5000', debug=True)
