# This builds the image locally for the current platform
# docker build --no-cache -t dynearthsol .

FROM ubuntu:22.04
LABEL maintainer="chaseshyu@gmail.com"

# Set environment variable to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install necessary software
RUN apt-get update && apt-get install -y \
      build-essential \
      python3 \
      python3-pip \
      vim \
      git \
      libboost-program-options-dev \ 
      && apt-get clean

# Add and enable the default user
ARG USER=human
ARG USER_ID=1000
RUN adduser --disabled-password --gecos '' --uid $USER_ID $USER
RUN adduser $USER sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Set $HOME
RUN chown -R $USER:$USER /home/$USER
USER $USER
ENV HOME /home/$USER
ENV USER $USER
WORKDIR $HOME

# pip install requirements
RUN pip install numpy

# git clone and build DynEarthSol
RUN git clone https://github.com/chaseshyu/DynEarthSol.git

# make dynearthsol2d
WORKDIR $HOME/DynEarthSol
RUN make ndims=2
WORKDIR $HOME
CMD ["bash"]
