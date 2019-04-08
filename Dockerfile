FROM nfcore/base
LABEL authors="harshil.patel@crick.ac.uk" \
      description="Docker image containing all requirements for the nfcore/atacseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-atacseq-1.0.0/bin:$PATH
