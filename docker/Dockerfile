FROM nanome/nanobabel

RUN apt install software-properties-common
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt update
RUN apt install -y python3.7
RUN apt install -y python3-pip

RUN pip3 install --upgrade pip
RUN python3.7 -m pip install nanome

COPY . /app
WORKDIR /app

CMD python3.7 run.py ${ARGS}