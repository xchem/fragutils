FROM python:3.13-slim

USER root
RUN apt-get --allow-releaseinfo-change update \
    && apt-get install -y \
        git \
        libfontconfig1 \
        libsm6 \
        libxrender1 \
        procps \
    && pip install rdkit==2025.3.6 \
    && git clone https://github.com/rdkit/mmpdb /usr/local/mmpdb \
    && pip install /usr/local/mmpdb

ADD . /usr/local/fragutils
RUN pip install /usr/local/fragutils
