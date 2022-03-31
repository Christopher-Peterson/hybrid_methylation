FROM rocker/tidyverse:4.1.2

# SeqMonk v1.48.1
# https://github.com/s-andrews/SeqMonk

# Setup timezone  
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Skipping the step for installing build environment; pretty sure it's already there

WORKDIR /usr/src

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/seqmonk_v1.48.1_linux64.tar.gz && \
  tar -xzf seqmonk_v1.48.1_linux64.tar.gz && \
  chmod a+x SeqMonk/seqmonk

ENV PATH=$PATH:/usr/src/SeqMonk/

# Figure out how to set the seqmonk scratch to the right directory