FROM tiangolo/uwsgi-nginx-flask:python3.8

# Copy requirements
COPY ./requirements.txt /var/www/requirements.txt
RUN pip install --upgrade pip \
    && pip install --default-timeout=1000 --no-cache-dir -r /var/www/requirements.txt

# Copy app code
COPY ./app /app
WORKDIR /app

# Run the Flask app
CMD [ "python", "run.py", "--host=0.0.0.0" ]