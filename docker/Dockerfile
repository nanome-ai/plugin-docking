FROM nanome/plugin-env

ENV ARGS=''
WORKDIR /app

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

# Set up conda environment for autodock4
RUN conda env create -f nanome_docking/autodock4/vina.yml

ENV SMINA_BINARY=nanome_docking/smina/smina_binary
RUN chmod +x  $SMINA_BINARY
CMD python run.py ${ARGS}
