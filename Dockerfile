FROM continuumio/miniconda3
MAINTAINER Patrick Hüther <patrick.huether@gmi.oeaw.ac.at>
LABEL authors="patrick.huether@gmi.oeaw.ac.at" \
    description="Container image containing all dependencies for the MethylScore pipeline"

#Install dependencies
RUN apt-get update 
RUN apt-get install -y build-essential procps libgsl-dev perl cpanminus
RUN apt-get clean -y

COPY environment.yaml /

#Install perl modules
RUN cpanm --notest Getopt::Long \
	Thread::Pool \
	Thread::Conveyor::Array \
	Thread::Conveyor::Throttled \
	Thread::Conveyor::Tied \
	Thread::Tie::Array \
	Config::Simple \
	File::Path \
	File::Copy \
	File::Tee \
	File::Which \
	Algorithm::Cluster \
	Data::Dumper \
	List::Util \
	Compress::BGZF::Reader

RUN conda env create -f /environment.yaml && conda clean -a
ENV PATH /opt/conda/envs/MethylScore/bin:$PATH

