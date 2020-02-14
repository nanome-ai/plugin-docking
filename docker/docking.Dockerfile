FROM continuumio/miniconda3

ENV ARGS=''

COPY . /app
WORKDIR /app

RUN pip install nanome
RUN conda install -c openbabel openbabel

CMD python run.py -a ${ARGS}
