FROM tiangolo/uwsgi-nginx-flask:python3.8
COPY ./requirements.txt /var/www/requirements.txt
RUN pip install --upgrade pip \
    && pip install --default-timeout=1000 --no-cache-dir -r /var/www/requirements.txt
COPY ./app/ /app
WORKDIR /app
CMD [ "python", "run.py", "--host=0.0.0.0"]