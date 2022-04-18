FROM ubuntu:20.04
# Runs Bedtools

# Deal w/ the time zone issue that can lag installs
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get -y upgrade && \
  apt-get install bedtools
