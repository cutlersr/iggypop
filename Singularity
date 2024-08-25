%post
    # Install necessary system packages for your application
    apt-get update && apt-get install -y \
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

    # Change to the /app directory
    mkdir -p /app
    cd /app

    # Copy all files to /app
    cp -r /context/* /app/

    # Update pip, setuptools, and wheel to the latest versions
    pip install --upgrade pip setuptools wheel

    # Install Python dependencies
    pip install -r /app/requirements.txt

    # Run setup script
    chmod +x /app/setup.sh
    /app/setup.sh

    # Ensure iggypop.py is executable
    chmod +x /app/iggypop.py

    # Create a symbolic link in /usr/local/bin
    ln -s /app/iggypop.py /usr/local/bin/iggypop
