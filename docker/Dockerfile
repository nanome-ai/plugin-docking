FROM nanome/plugin-env

ENV ARGS=''
WORKDIR /app

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt

# adfr-suite functions required for autodock4. Requires separate python 2.7 environment.
COPY adfr-suite.yml .
RUN conda env create --file adfr-suite.yml

COPY . .

# Openbabel also required for autodock4
RUN conda install -c conda-forge openbabel

# Ensure executables are executable.
ENV SMINA_BINARY=nanome_docking/smina/smina_binary
ENV VINA_BINARY=nanome_docking/autodock4/vina_1.2.2_linux_x86_64
RUN chmod +x  $SMINA_BINARY
RUN chmod +x  $VINA_BINARY
CMD python run.py ${ARGS}
