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
        git
        swig
        libssl-dev
        libffi-dev
        libxml2-dev
        libxslt1-dev
        zlib1g-dev
        r-base
        libhdf5-dev
        wget
        ncbi-blast+
        bowtie
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

# Ensure mfeprimer is executable
chmod +x data/mfeprimer-3.3.1-linux-386

# Determine the directory of the current script
SCRIPT_DIR=$(dirname "$(realpath "$0")")

# Ensure iggypop.py is executable
chmod +x "$SCRIPT_DIR/iggypop.py"

# Create ~/.local/bin if it doesn't exist and ensure it's in the PATH
mkdir -p ~/.local/bin
if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
    echo 'export PATH=$PATH:$HOME/.local/bin' >> ~/.bashrc
    export PATH=$PATH:$HOME/.local/bin
fi

# Create a symbolic link in ~/.local/bin
if [ ! -e ~/.local/bin/iggypop ]; then
    ln -s "$SCRIPT_DIR/iggypop.py" ~/.local/bin/iggypop
fi

echo "Setup complete. Please restart your terminal or run 'source ~/.bashrc' to update your PATH if necessary."
