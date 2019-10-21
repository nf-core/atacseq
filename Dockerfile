FROM nfcore/base:1.7
LABEL authors="Harshil Patel" \
      description="Docker image containing all requirements for nf-core/atacseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-atacseq-1.0.1dev > nf-core-atacseq-1.0.1dev.yml
ENV PATH /opt/conda/envs/nf-core-atacseq-1.0.1dev/bin:$PATH
