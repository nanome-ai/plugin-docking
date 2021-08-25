FROM nanome/plugin-env

ENV ARGS=''
WORKDIR /app

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

ENV SMINA_EXEC=smina/smina
RUN chmod +x  $SMINA_EXEC
CMD python run.py ${ARGS}
