from flask import Flask, send_from_directory
import os
import redis

app = Flask(__name__)
app.config['SECRET_KEY'] = '5791628bb0b13ce0c676dfde280ba245'

# Redis configuration
REDIS_HOST = os.environ.get('REDIS_HOST', 'localhost')
REDIS_PORT = int(os.environ.get('REDIS_PORT', 6379))
REDIS_DB = int(os.environ.get('REDIS_DB', 0))

try:
    redis_client = redis.Redis(
        host=REDIS_HOST,
        port=REDIS_PORT,
        db=REDIS_DB,
        decode_responses=True  # This helps work with strings more easily
    )
    # Optionally test the connection
    redis_client.ping()
except redis.ConnectionError:
    print("Failed to connect to Redis. Ensure Redis server is running.")
    redis_client = None

@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'), 'favicon.ico', mimetype='image/vnd.microsoft.icon')

@app.route('/static/data/<filename>')
def serve_mapping_files(filename):
    return send_from_directory(os.path.join(app.root_path, 'static', 'data'), filename)

# LAUNCH the app
from routes import *
