# This builds the image locally for the current platform
# docker pull ubuntu:22.04
# docker build -t dynearthsol .
# docker run --rm -it dynearthsol

# Use the Ubuntu 22.04 release as the base image
ARG BASE_IMAGE=ubuntu:22.04
FROM ${BASE_IMAGE}
# Set the maintainer
LABEL maintainer="chaseshyu@gmail.com"

# Set environment variable to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

ARG TZ=UTC
ARG GCCVERSION=gcc-11
# Update package list and install necessary software
RUN apt-get update && apt-get install -y \
      build-essential \
      libboost-program-options-dev \
      python3 \
      python3-pip \
      vim \
    && if [ "$GCCVERSION" = "gcc-8" ]; then \
      # Install gcc-8 and g++-8
      apt-get install -y software-properties-common \
      && add-apt-repository ppa:ubuntu-toolchain-r/test \
      && apt-get update \
      && apt-get install -y gcc-8 g++-8 \
      && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 60 \
      && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-8 60; \
    elif [ "$GCCVERSION" = "gcc-11" ]; then \
      apt-get install -y tzdata; \
    fi \
    && apt-get clean \
    # set time zone of system
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# Add and enable the default user
ARG USER=human
ARG USER_ID=1000
RUN adduser --disabled-password --gecos '' --uid $USER_ID $USER \
    && adduser $USER sudo \
    && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers \
    && chown -R $USER:$USER /home/$USER

# Set $HOME
USER $USER
ENV HOME /home/$USER
ENV USER $USER

# pip install requirements
RUN pip install numpy scipy

# Clone the repository
COPY --chown=$USER:$USER . $HOME/DynEarthSol

# make dynearthsol2d
WORKDIR $HOME/DynEarthSol
ARG NDIMS=2
RUN make cleanall && make ndims=$NDIMS -j4

# Default use 8 threads
ENV OMP_NUM_THREADS=8

CMD ["bash"]