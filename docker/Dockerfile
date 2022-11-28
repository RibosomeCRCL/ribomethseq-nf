FROM continuumio/miniconda3:4.12.0

COPY conda.yml .
COPY Rscripts/ ./Rscripts/
COPY ref_human/ ./ref_human/
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

RUN apt-get install -y procps

