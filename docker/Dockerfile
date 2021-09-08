FROM nanome/plugin-env

ENV ARGS=''
WORKDIR /app

ARG CACHEBUST
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

# Set up conda environment for adfr-suite functions for vina prep steps that require python 2.7
RUN conda env create --file /app/nanome_docking/autodock4/adfr-suite.yml

RUN conda install -c conda-forge openbabel

# Ensure executables are executable.
ENV SMINA_BINARY=nanome_docking/smina/smina_binary
ENV VINA_BINARY=nanome_docking/autodock4/vina_1.2.2_linux_x86_64
RUN chmod +x  $SMINA_BINARY
RUN chmod +x  $VINA_BINARY
CMD python run.py ${ARGS}
