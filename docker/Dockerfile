FROM nanome/plugin-env

ENV ARGS=''
WORKDIR /app

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

ENV SMINA_BINARY=smina/smina_binary
RUN chmod +x  $SMINA_BINARY
CMD python run.py ${ARGS}
