FROM nanome/plugin-env

ENV ARGS=''
WORKDIR /app

ARG CACHEBUST

# Set Arg to anything to avoid installing Autodock4 dependencies.
ARG EXCLUDE_AUTODOCK4

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY adfr-suite.yml .

RUN if [ ! -z "$EXCLUDE_AUTODOCK4" ]; \
    then \
    echo "Installing Autodock4 requirements from Conda." \
    # adfr-suite functions required for autodock4. Requires separate python 2.7 environment.
    && conda env create --file adfr-suite.yml \
    # Openbabel also required for autodock4
    && conda install -c conda-forge openbabel; \
    else \
    echo "Not installing Autodock4 requirements"; \
    fi

COPY . .

# Ensure executables are executable.
ENV SMINA_BINARY=nanome_docking/smina/smina_binary
ENV VINA_BINARY=nanome_docking/autodock4/vina_1.2.2_linux_x86_64
RUN chmod +x $SMINA_BINARY
RUN chmod +x $VINA_BINARY
CMD python run.py ${ARGS}
