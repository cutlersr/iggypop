# Singularity definition file

BootStrap: docker
From: python:3.9-slim

%labels
    Maintainer "cutler@ucr.edu"

%post
    # Update the package list and install necessary system packages
    apt-get update && apt-get install -y \
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

    # Update pip, setuptools, and wheel to the latest versions and install Python packages
    pip install --upgrade pip setuptools wheel
    pip install -r /app/requirements.txt

    # Ensure setup script is executable and run it
    chmod +x /app/setup.sh && ./app/setup.sh

    # Make the Python script executable
    chmod +x /app/iggypop.py

%files
    # Copy all files to /app
    . /app

%environment
    # Set the working directory
    export SINGULARITY_WORKDIR=/app

%runscript
    exec python3 /app/iggypop.py
