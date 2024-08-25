# Use an official Python runtime as a parent image with amd64 architecture
FROM python:3.9-slim

# Set a maintainer label
LABEL maintainer="cutler@ucr.edu"

# Install necessary system packages for your application
RUN apt-get update && apt-get install -y \
    git \
    swig \
    libssl-dev \
    libffi-dev \
    libxml2-dev \
    libxslt1-dev \
    zlib1g-dev \
    r-base \
    libhdf5-dev \
    wget \
    ncbi-blast+ \
    bowtie && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /app
COPY . /app

# Update pip, setuptools, and wheel to the latest versions
RUN pip install --upgrade pip setuptools wheel && \
    pip install -r requirements.txt

# Run setup script
RUN chmod +x setup.sh && ./setup.sh

# Create a symbolic link to the iggypop.py script in /usr/local/bin
RUN ln -s /app/iggypop.py /usr/local/bin/iggypop

# Ensure the script has executable permissions
RUN chmod +x /app/iggypop.py

# Set the default command to run the script
CMD ["iggypop"]
