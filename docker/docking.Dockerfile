FROM python:3.7

ENV PLUGIN_SERVER=plugins.nanome.ai

COPY . /app
WORKDIR /app

RUN pip install nanome
RUN chmod +x nanome_docking/smina

CMD python -m nanome_docking.Docking smina -a ${PLUGIN_SERVER}