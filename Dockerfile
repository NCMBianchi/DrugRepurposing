FROM tiangolo/uwsgi-nginx-flask:python3.8

# INSTALLATION
RUN apt-get update && apt-get install -y \
    docker.io \
    docker-compose \
    && rm -rf /var/lib/apt/lists/*

# PARAMETERS AND ENTRYPOINT (1)
COPY Docker_parameters.txt /app/Docker_parameters.txt
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# PYTHON
COPY requirements.txt /var/www/requirements.txt
RUN pip install --upgrade pip \
    && pip install --default-timeout=1000 --no-cache-dir -r /var/www/requirements.txt \
    && pip install docker

# FILES
COPY ./app/ /app
WORKDIR /app

# ENTRYPOINT (2)
ENTRYPOINT ["/entrypoint.sh"]
CMD [ "python", "run.py", "--host=0.0.0.0"]