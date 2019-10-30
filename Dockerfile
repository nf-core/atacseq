FROM nfcore/base:1.7
LABEL authors="Harshil Patel" \
      description="Docker image containing all requirements for nf-core/atacseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-atacseq-1.1.0 > nf-core-atacseq-1.1.0.yml
ENV PATH /opt/conda/envs/nf-core-atacseq-1.1.0/bin:$PATH
