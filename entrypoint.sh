#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "Starting DrugApp main service..."

# Print Python version and environment info
python --version
echo "Working directory: $(pwd)"
echo "Contents of /app:"
ls -la /app

# Start the Flask application with more verbose output
exec uwsgi --http "0.0.0.0:5000" --module run:app --need-app --python-autoreload=1