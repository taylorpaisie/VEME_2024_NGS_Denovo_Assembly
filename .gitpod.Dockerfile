FROM gitpod/workspace-full

USER gitpod

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Initialize Conda
RUN $HOME/miniconda/bin/conda init bash

# Add the following lines to your .gitpod.Dockerfile file:
FROM gitpod/workspace-full-vnc:latest
USER root
RUN apt-get update && apt-get install -y x11-apps
USER gitpod
