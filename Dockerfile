FROM python:3.8

ADD main.py .

ADD statistics.py .

ADD genePop5Ix5L .

ADD OneSamp .

ADD ./scripts/rScript.r .

RUN pip install numpy

ENTRYPOINT ["python", "./main.py"]
