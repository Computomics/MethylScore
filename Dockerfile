FROM debian:bullseye-slim AS BUILDER
LABEL authors="patrick.huether@gmi.oeaw.ac.at" \
    description="Container image containing all dependencies for the MethylScore pipeline"

#Install dependencies 
RUN apt-get update && apt-get install --no-install-recommends -y git ca-certificates build-essential libgsl-dev perl cpanminus

#Install perl modules
RUN cpanm Algorithm::Cluster Compress::BGZF::Reader

#Compile hmm_mrs
RUN git clone --depth 1 --branch master https://github.com/Computomics/MethylScore

WORKDIR /MethylScore/src/hmm_mrs/src

RUN cd smithlab_cpp \
    && make \
    && cd ../analysis \
    && SMITHLAB_CPP=../smithlab_cpp make \
    && chmod +x hmm_mrs

FROM mambaorg/micromamba:0.19.1

COPY --from=BUILDER /usr/local /usr/local
COPY --from=BUILDER /MethylScore/src/hmm_mrs/src/analysis/hmm_mrs /MethylScore/bin/betabin_model /usr/local/bin/
COPY --chown=micromamba:micromamba environment.yaml /tmp/environment.yaml

USER root
RUN apt-get update \
	&& apt-get install --no-install-recommends -y \
	procps perl libconfig-simple-perl libfile-which-perl libfile-tee-perl libthread-conveyor-perl libthread-pool-perl \
	&& apt-get clean -y \
	&& rm -rf /var/lib/{apt,dpkg,cache,log}
USER micromamba

RUN micromamba install -y -n base -f /tmp/environment.yaml && micromamba clean -a
