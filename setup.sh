#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Function to install R packages if not already installed
install_r_package_if_needed() {
    Rscript -e "if (!requireNamespace('$1', quietly = TRUE)) install.packages('$1', repos='https://cloud.r-project.org')"
}

# Function to check and install packages using Conda
check_and_install_conda() {
    local PACKAGE=$1
    if ! conda list | grep -q "^$PACKAGE "; then
        echo "$PACKAGE is not installed. Installing via Conda..."
        conda install -y -c bioconda -c conda-forge "$PACKAGE"
    else
        echo "$PACKAGE is already installed."
    fi
}

# Determine if running inside Docker/Singularity 
if [ -f /.dockerenv ]; then
    echo "Running inside a Docker container. Skipping system package installation."
elif grep -qE 'singularity' /proc/self/cgroup; then
    echo "Running inside a Singularity container. Skipping system package installation."
else
    echo "Installing system dependencies via Conda..."

    CONDA_PACKAGES=(
        bowtie        
        git
        bzip2          
        expat          
        libffi          
        gdbm          
        krb5            
        hdf5            
        keyutils        
        xz             
        ncurses         
        libnsl         
        readline       
        sqlite         
        openssl         
        libtirpc        
        libuuid         
        libxml2         
        libxslt        
        blast          
        r-base
        swig
        tzdata
        wget
        zlib            
        libxcrypt       
        libdb           
        )


    for PACKAGE in "${CONDA_PACKAGES[@]}"; do
        check_and_install_conda "$PACKAGE"
    done
fi

# Install Python packages
echo "Installing Python packages..."
pip install --no-cache-dir -r requirements.txt

# Install R packages
echo "Installing R packages..."
Rscript install_packages.R

SCRIPT_DIR=$(dirname "$(realpath "$0")")
# Ensure iggypop.py is executable
chmod +x "$SCRIPT_DIR/iggypop.py"

# Download the donor_200.h5 and acceptor_200.h5 files to the data/ folder
echo "Downloading donor_200.h5 and acceptor_200.h5 to the data/ folder..."

DATA_DIR="$SCRIPT_DIR/data"
mkdir -p "$DATA_DIR"

wget -O "$DATA_DIR/donor_200.h5" "https://git.unistra.fr/nscalzitti/spliceator/-/raw/master/Models/donor_200.h5"
wget -O "$DATA_DIR/acceptor_200.h5" "https://git.unistra.fr/nscalzitti/spliceator/-/raw/master/Models/acceptor_200.h5"

echo "Files downloaded to the data/ folder."

echo "Setup complete. Please restart your terminal or run 'source ~/.bashrc' to update your PATH if necessary."
