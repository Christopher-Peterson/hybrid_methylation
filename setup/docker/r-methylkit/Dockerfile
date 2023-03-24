FROM crpeters/r-plotting:4.2.1

RUN install2.r -e BiocManager
# Install a pre-req
RUN apt-get update && apt-get install -y liblzma-dev libbz2-dev libcurl4-openssl-dev

RUN /usr/local/lib/R/site-library/littler/examples/installBioc.r Rhtslib

RUN /usr/local/lib/R/site-library/littler/examples/installBioc.r methylKit DSS