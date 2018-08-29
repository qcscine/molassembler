# Inspired by https://github.com/orian/cppenv/blob/master/Dockerfile
FROM ubuntu:bionic

LABEL description="A linux C++ build environment on ubuntu 18.04 with GCC 7."

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  ca-certificates \
  cmake \
  gcc \
  git \
  libeigen3-dev \
  liblapack-dev \
  libopenblas-dev \
  make \
  wget

# Install boost
COPY install_boost.sh /install_boost.sh
RUN bash install_boost.sh
ENV BOOST_ROOT=/soft/boost

# Install dlib
COPY install_dlib.sh /install_dlib.sh
RUN bash install_dlib.sh
ENV dlib_DIR=/soft/dlib/lib/cmake/dlib

# Install scine
COPY scine.tar /scine.tar
COPY install_scine.sh /install_scine.sh
RUN bash install_scine.sh
ENV scine_DIR=/soft/scine/lib/cmake/scine

# remove some unneeded weight
COPY cleanup.sh /cleanup.sh
RUN bash cleanup.sh

# setup working directory
RUN mkdir /devel
WORKDIR /devel

CMD ["/bin/bash"]
