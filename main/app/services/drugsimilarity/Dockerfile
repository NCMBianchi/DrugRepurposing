# Drug Similarity Service Dockerfile
FROM python:3.8-slim-buster

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY app/ .

# Environment variables for configuration
ENV INPUT_SEED=""
ENV DATE=""
ENV BASE_DIRECTORY=""
ENV MIN_SIMILARITY=""
ENV NODES=""

# Default command
CMD ["python", "drugsimilarity_service.py"]