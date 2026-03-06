FROM python:3.11.15-slim

USER root
RUN apt-get --allow-releaseinfo-change update \
    && apt-get install -y \
        git \
        libfontconfig1 \
        libsm6 \
        libxrender1 \
        procps \
    && pip install rdkit==2023.3.2 \
    && git clone https://github.com/rdkit/mmpdb /usr/local/mmpdb \
    && pip install /usr/local/mmpdb

ADD . /usr/local/frag
RUN pip install /usr/local/frag
