FROM continuumio/miniconda3:4.12.0

LABEL authors="Théo COMBE" \
      description="Docker image for ribomethseq-nf workflow"

COPY conda.yml .

RUN    conda env update -n root -f conda.yml \
    && conda clean -a

# Nextflow needs this
RUN apt-get install -y procps

# Dump the details of the installed packages to a file for posterity
RUN conda env export > ribomethseq-nf.yml
