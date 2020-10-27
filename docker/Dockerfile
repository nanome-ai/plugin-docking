FROM nanome/plugin-env

ENV ARGS=''

ARG CACHEBUST
RUN pip3 install nanome

COPY . /app
WORKDIR /app

CMD python3 run.py ${ARGS}
