FROM nfcore/base
LABEL authors="Harshil Patel" \
      description="Docker image containing all requirements for nf-core/atacseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-atacseq-1.0.1dev/bin:$PATH
