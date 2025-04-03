# Use an official Python runtime as a parent image with amd64 architecture
FROM python:3.9-slim

# Set a maintainer label
LABEL maintainer="cutler@ucr.edu"

# Install necessary system packages for your application
RUN apt-get update && apt-get install -y \
    bowtie \
    git \
    libbz2-1.0 \
    libc6 \
    libcom-err2 \
    libcrypt1 \
    libdb5.3 \
    libexpat1 \
    libffi-dev \
    libffi8 \
    libgdbm6 \
    libgssapi-krb5-2 \
    libhdf5-dev \
    libk5crypto3 \
    libkeyutils1 \
    libkrb5-3 \
    libkrb5support0 \
    liblzma5 \
    libncursesw6 \
    libnsl2 \
    libreadline8 \
    libsqlite3-0 \
    libssl-dev \
    libssl3 \
    libtinfo6 \
    libtirpc3 \
    libuuid1 \
    libxml2-dev \
    libxslt1-dev \
    ncbi-blast+ \
    netbase \
    r-base \
    swig \
    tzdata \
    wget \
    zlib1g \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda to provide the 'conda' command
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add Conda to PATH so that all subsequent RUN commands have access to it
ENV PATH=/opt/conda/bin:$PATH

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Update pip, setuptools, and wheel to the latest versions and install Python dependencies
RUN pip install --upgrade pip setuptools wheel && \
    pip install -r requirements.txt

# Run the setup script (which uses conda); note that it should now find conda in PATH
RUN chmod +x setup.sh && ./setup.sh

# Ensure the main application script has executable permissions
RUN chmod +x /app/iggypop.py

# Set the container's entry point to bash
ENTRYPOINT ["/bin/bash"]
