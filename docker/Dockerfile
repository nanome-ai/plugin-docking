ARG ALGORITHM=smina

# Multiple build options based on ALGORITHM build-arg
# Avoids installing uneecessary dependencies and long build times

# This part runs for every algorithm install
FROM nanome/plugin-env AS base
ENV ARGS=''
WORKDIR /app
COPY requirements.txt .
ENV ALGORITHM=smina
ARG CACHEBUST
RUN pip install -r requirements.txt

COPY . .

# Section of build that runs when smina selected (Default)
FROM base AS smina-build
RUN echo "ALGORITHM=smina"
ENV SMINA_BINARY=plugin/smina/smina_binary
RUN chmod +x $SMINA_BINARY

# Section of build that runs when autodock4 selected
FROM base AS autodock4-build
RUN echo "ALGORITHM=autodock4"
ENV ALGORITHM=autodock4
ENV VINA_BINARY=plugin/autodock4/vina_1.2.2_linux_x86_64
COPY adfr-suite.yml .
RUN conda env create --file adfr-suite.yml
RUN conda install -c conda-forge openbabel
RUN chmod +x $VINA_BINARY

# Determine which branch to build based on build-arg.
FROM ${ALGORITHM}-build
CMD python run.py --algorithm ${ALGORITHM} ${ARGS}

