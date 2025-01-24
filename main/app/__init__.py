from flask import Flask, send_from_directory
import os

app = Flask(__name__)
app.config['SECRET_KEY'] = '5791628bb0b13ce0c676dfde280ba245'

@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'), 'favicon.ico', mimetype='image/vnd.microsoft.icon')

@app.route('/static/data/<filename>')
def serve_mapping_files(filename):
    return send_from_directory(os.path.join(app.root_path, 'static', 'data'), filename)

# LAUNCH the app
from routes import *