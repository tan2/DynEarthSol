# This builds the image locally for the current platform
# docker pull ubuntu:22.04
# docker build -t dynearthsol .
# docker run --rm -it dynearthsol

FROM ubuntu:22.04
LABEL maintainer="chaseshyu@gmail.com"

# Set environment variable to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
# Default use 8 threads
ENV OMP_NUM_THREADS=8

# Update package list and install necessary software
RUN apt-get update && apt-get install -y \
      build-essential libboost-program-options-dev \ 
      python3 python3-pip vim git \
      && apt-get clean

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
WORKDIR $HOME

# pip install requirements
RUN pip install numpy scipy

# Clone the repository
COPY --chown=$USER:$USER . $HOME/DynEarthSol

# make dynearthsol2d
WORKDIR $HOME/DynEarthSol
RUN make cleanall && make ndims=2 -j4

CMD ["bash"]