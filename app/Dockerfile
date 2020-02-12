FROM python:3.8-buster

MAINTAINER Zi Han Ooi "bt19623@qmul.ac.uk"

RUN apt update && apt install build-essential

COPY ./requirements.txt /app/requirements.txt

WORKDIR /app

RUN pip install -r requirements.txt

COPY . /app

EXPOSE 5000

CMD [ "python3", "app.py" ]
