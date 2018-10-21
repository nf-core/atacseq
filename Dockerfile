FROM nfcore/base

MAINTAINER Harshil Patel <harshil.patel@crick.ac.uk>
LABEL authors="harshil.patel@crick.ac.uk" \
    description="Docker image containing all requirements for the nfcore/atacseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-atacseq-1.0dev/bin:$PATH
