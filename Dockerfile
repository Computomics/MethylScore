FROM continuumio/miniconda3
MAINTAINER Patrick Hüther <patrick.huether@gmi.oeaw.ac.at>
LABEL authors="patrick.huether@gmi.oeaw.ac.at" \
    description="Container image containing all dependencies for the MethylScore pipeline"

COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/MethylScore/bin:$PATH
