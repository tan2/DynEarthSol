#!/bin/bash

# Set the number of dimensions
NDIMS=2 # or 3
# Set the GCC version
GCCVERSION=gcc-11 # or gcc-8
# Set the timezone of the host machine
HOST_TZ=$(readlink /etc/localtime | sed 's|.*/zoneinfo/||')

# based on gcc version to choose the image
if [ "$GCCVERSION" = "gcc-8" ]; then
    BASE_IMAGE=ubuntu:20.04
elif [ "$GCCVERSION" = "gcc-11" ]; then
    BASE_IMAGE=ubuntu:22.04
else
    echo "GCC version not supported"
    exit 1
fi

# Build the docker image
docker build -t dynearthsol/$GCCVERSION . \
    --build-arg GCCVERSION=$GCCVERSION \
    --build-arg BASE_IMAGE=$BASE_IMAGE \
    --build-arg NDIMS=$NDIMS \
    --build-arg TZ=$HOST_TZ \

echo "To run the container execute:"
echo "$ docker run -it --rm dynearthsol/$GCCVERSION"