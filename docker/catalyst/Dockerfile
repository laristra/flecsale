
# Download base image ubuntu 16.04
FROM ubuntu:devel

ENV http_proxy=http://proxyout.lanl.gov:8080
ENV https_proxy=http://proxyout.lanl.gov:8080

# Author
MAINTAINER Pascal Grosset

# Define default command.
CMD ["bash"]

# Update and install packages needed
RUN apt-get update
RUN apt-get -y install git wget cmake openssh-server build-essential python2.7 python-dev libopenmpi-dev openmpi-bin freeglut3-dev qt4-dev-tools

# Expose the ParaView Catalyst Port
EXPOSE 22222


# Set the default working directory for the container
WORKDIR /home


# Go to home and download Paraview and Catalyst flecsale
RUN wget https://www.paraview.org/files/v5.2/ParaView-v5.2.0.tar.gz 
RUN tar -zxvf ParaView-v5.2.0.tar.gz

RUN git clone -b bugfix/vtk --recursive https://github.com/laristra/flecsale.git
RUN ls


# Build ParaView 5.2
WORKDIR /home/ParaView-v5.2.0
RUN mkdir build
WORKDIR /home/ParaView-v5.2.0/build
RUN ls
RUN cmake .. -DPARAVIEW_USE_MPI:BOOL=ON -DPARAVIEW_ENABLE_PYTHON:BOOL=ON -DPARAVIEW_ENABLE_CATALYST:BOOL=ON
RUN make -j16
ENV ParaView_DIR="/home/ParaView-v5.2.0/build"



# Build Flecsale Catalyst
WORKDIR /home/flecsale
RUN mkdir build
WORKDIR /home/flecsale/build
RUN cmake .. -DENABLE_CATALYST:BOOL=ON -DCMAKE_DISABLE_FIND_PACKAGE_VTK=on
RUN make 
