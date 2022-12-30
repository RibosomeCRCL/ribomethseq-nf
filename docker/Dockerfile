###
### Building all required tools
###
FROM debian:bullseye-20221219-slim as builder

RUN mkdir -p /home/rms/bin
WORKDIR /home/rms/src

#
# R
#
RUN    apt-get update                                \
    && apt-get install -y wget gcc g++ gfortran make \
        cmake pkg-config xorg-dev xvfb               \
        libjpeg-dev libpng-dev libtiff-dev           \
        libc6-dev libbz2-dev                         \
        libcurl4-openssl-dev libgomp1 liblzma-dev    \
        libpcre++-dev libpcre2-dev libquadmath0      \
        libreadline-dev zlib1g-dev libcairo2-dev     \
        libxml2-dev

RUN    wget https://cloud.r-project.org/src/base/R-4/R-4.2.2.tar.gz \
    && tar xzf R-4.2.2.tar.gz

RUN    cd R-4.2.2 \
    && ./configure --prefix=/usr/R --without-aqua --disable-nls --disable-java \
    && make all && make install

ENV PATH="/usr/R/bin:${PATH}"

# factoextra imports tidyr, reshape2
# tidyr imports dplyr, tibble
RUN R -e "install.packages('ade4',       dependencies=c('Depends', 'Imports', 'LinkingTo'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('factoextra', dependencies=c('Depends', 'Imports', 'LinkingTo'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('rmarkdown',  dependencies=c('Depends', 'Imports', 'LinkingTo'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2',    dependencies=c('Depends', 'Imports', 'LinkingTo'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('pheatmap',   dependencies=c('Depends', 'Imports', 'LinkingTo'), repos='http://cran.rstudio.com/')"


# more requirements
RUN    apt-get update                                           \
    && apt-get install -y autoconf automake perl libdeflate-dev \
                          unzip zip default-jre python3 pip

#
# htslib / samtools
#

RUN    wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 \
    && tar xjf htslib-1.16.tar.bz2

RUN    cd htslib-1.16 && autoreconf -i                 \
    && ./configure --disable-plugins --disable-libcurl \
    && make -j4 libhts.a bgzip htsfile tabix           \
    && mv bgzip htsfile tabix /home/rms/bin

RUN    wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
    && tar xjf samtools-1.16.1.tar.bz2

RUN    cd samtools-1.16.1 && autoheader && autoconf -Wno-syntax                     \
    && ./configure --without-curses --without-ncursesw --with-htslib=../htslib-1.16 \
    && make && mv samtools /home/rms/bin

#
# bedtools
#
RUN    wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
    && chmod +x bedtools.static.binary \
    && mv bedtools.static.binary /home/rms/bin/bedtools

#
# bowtie2
#
RUN    wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.0/bowtie2-2.5.0-linux-x86_64.zip \
    && unzip bowtie2-2.5.0-linux-x86_64.zip \
    && cd bowtie2-2.5.0-linux-x86_64 \
    && mv bowtie2 bowtie2-align-s bowtie2-align-l       \
          bowtie2-build bowtie2-build-s bowtie2-build-l \
          bowtie2-inspect bowtie2-inspect-s bowtie2-inspect-l /home/rms/bin


RUN    cd /home/rms/bin \
    && strip bowtie2-align-s bowtie2-align-l bowtie2-build-s bowtie2-build-l \
             bowtie2-inspect-s bowtie2-inspect-l bedtools samtools tabix bgzip htsfile

#
# Trimmomatic
#
RUN    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip \
    && unzip Trimmomatic-0.39.zip \
    && mkdir -p /home/rms/Trimmomatic \
    && mv Trimmomatic-0.39/trimmomatic-0.39.jar /home/rms/Trimmomatic

COPY trimmomatic /home/rms/bin

#
# FastQC
#
RUN    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip \
    && chmod +x FastQC/fastqc && mv FastQC /home/rms \
    && cd /home/rms/bin && ln -s ../FastQC/fastqc

#
# MultiQC
#
RUN pip install multiqc

#
# pandoc
#
RUN    wget https://github.com/jgm/pandoc/releases/download/2.19.2/pandoc-2.19.2-1-amd64.deb \
    && dpkg -i pandoc-2.19.2-1-amd64.deb \
    && mv /usr/bin/pandoc /usr/local/bin

ENV PATH="/home/rms/bin:${PATH}"

WORKDIR /home/rms

CMD ["/bin/bash"]
