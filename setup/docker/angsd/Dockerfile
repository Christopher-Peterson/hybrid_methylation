FROM crpeters/trim_galore:0.6.7
# ANGSD 0.940-stable
# Download & extract the package
RUN wget https://github.com/ANGSD/angsd/releases/download/0.940/angsd0.940.tar.gz && \
  tar xf angsd0.940.tar.gz
# Install a pre-req
RUN apt-get update && apt-get install -y liblzma-dev libbz2-dev libcurl4-openssl-dev

RUN cd htslib && make

RUN cd angsd && \
  make HTSSRC=../htslib

RUN mkdir /angsd-files && mv angsd /angsd-files/ && mv htslib /angsd-files

ENV PATH=$PATH:/angsd-files/angsd
