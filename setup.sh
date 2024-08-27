#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Function to install R packages if not already installed
install_r_package_if_needed() {
    Rscript -e "if (!requireNamespace('$1', quietly = TRUE)) install.packages('$1', repos='https://cloud.r-project.org')"
}

# Function to check if a package is installed, and install it if it is not
check_and_install() {
    local PACKAGE=$1
    if ! dpkg -l | grep -q "^ii  $PACKAGE "; then
        echo "$PACKAGE is not installed. Installing..."
        apt-get install -y "$PACKAGE"
    else
        echo "$PACKAGE is already installed."
    fi
}

# Determine if running inside Docker or Singularity
if [ -f /.dockerenv ]; then
    echo "Running inside a Docker container. Skipping system package installation."
elif grep -qE 'singularity' /proc/self/cgroup; then
    echo "Running inside a Singularity container. Skipping system package installation."
else
    # Ensure necessary system packages are installed
    echo "Checking and installing necessary system packages..."
    REQUIRED_PACKAGES=(
    bowtie 
    git 
    libbz2-1.0 
    libc6 
    libcom-err2 
    libcrypt1 
    libdb5.3 
    libexpat1 
    libffi-dev 
    libffi8 
    libgdbm6 
    libgssapi-krb5-2 
    libhdf5-dev 
    libk5crypto3 
    libkeyutils1 
    libkrb5-3 
    libkrb5support0 
    liblzma5 
    libncursesw6 
    libnsl2 
    libreadline8 
    libsqlite3-0 
    libssl-dev 
    libssl3 
    libtinfo6 
    libtirpc3 
    libuuid1 
    libxml2-dev 
    libxslt1-dev 
    ncbi-blast+ 
    netbase 
    r-base 
    swig 
    tzdata 
    wget 
    zlib1g 
    zlib1g-dev
    )

    for PACKAGE in "${REQUIRED_PACKAGES[@]}"; do
        check_and_install "$PACKAGE"
    done
fi

# Install Python packages
echo "Installing Python packages..."
pip install --no-cache-dir -r requirements.txt

# Install R packages
echo "Installing R packages..."
Rscript install_packages.R

# Determine the directory of the current script
SCRIPT_DIR=$(dirname "$(realpath "$0")")

# Ensure iggypop.py is executable
chmod +x "$SCRIPT_DIR/iggypop.py"

echo "Setup complete. Please restart your terminal or run 'source ~/.bashrc' to update your PATH if necessary."
