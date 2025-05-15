import os
from flask import Flask, render_template_string, request, redirect, url_for, flash
from werkzeug.utils import secure_filename

# Configuration
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'txt', 'csv', 'pdf', 'png', 'jpg', 'jpeg', 'gif'}

app = Flask(__name__)
app.secret_key = 'supersecretkey'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

# Example pipelines
PIPELINES = ['RNA-Seq', 'Variant Calling', 'Metagenomics']

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        pipeline = request.form.get('pipeline')
        file = request.files.get('file')
        if not pipeline:
            flash('Please select a pipeline.')
            return redirect(request.url)
        if not file or file.filename == '':
            flash('No file selected.')
            return redirect(request.url)
        if allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            flash(f"File '{filename}' uploaded for pipeline '{pipeline}'.")
            return redirect(request.url)
        else:
            flash('File type not allowed.')
            return redirect(request.url)
    # Simple HTML template for demonstration
    return render_template_string('''
        <!doctype html>
        <title>Bioinformatics pipelines</title>
        <h1>Bioinformatics pipelines</h1>
        {% with messages = get_flashed_messages() %}
          {% if messages %}
            <ul>
            {% for message in messages %}
              <li>{{ message }}</li>
            {% endfor %}
            </ul>
          {% endif %}
        {% endwith %}
        <form method=post enctype=multipart/form-data>
          <label for="pipeline">Choose pipeline:</label>
          <select name="pipeline" id="pipeline" required>
            <option value="" disabled selected>Select pipeline</option>
            {% for p in pipelines %}
              <option value="{{ p }}">{{ p }}</option>
            {% endfor %}
          </select><br><br>
          <input type=file name=file required>
          <input type=submit value=Upload>
        </form>
    ''', pipelines=PIPELINES)

if __name__ == '__main__':
    app.run(debug=True)
