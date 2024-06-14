FROM gitpod/workspace-full

USER gitpod

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Initialize Conda
RUN $HOME/miniconda/bin/conda init bash


# Use a base image with VNC and X11 tools
FROM gitpod/workspace-full-vnc:latest

# Install Miniconda and X11 tools
USER root
RUN apt-get update && apt-get install -y x11-apps

# Switch back to the gitpod user
USER gitpod

# Add any other dependencies or configurations here

