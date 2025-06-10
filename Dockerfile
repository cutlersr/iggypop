FROM python:3.9-slim
LABEL maintainer="cutler@ucr.edu"

# ---------- system deps ----------
RUN apt-get update && apt-get install -y \
      build-essential git wget swig tzdata r-base \
      libbz2-1.0 libexpat1 libffi-dev libgdbm6 \
      libxml2-dev libxslt1-dev libssl-dev zlib1g-dev libtirpc-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ---------- miniconda ----------
RUN wget -O /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh \
 && bash /tmp/miniconda.sh -b -p /opt/conda && rm /tmp/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# ---------- create env ----------
RUN conda create -y -n iggypop python=3.9 && conda clean -afy
ENV PATH=/opt/conda/envs/iggypop/bin:$PATH

# ---------- install packages in *that* env ----------
RUN conda install -y -n iggypop -c conda-forge -c bioconda \
      pandas numpy scipy cython bowtie \
      git bzip2 expat libffi gdbm krb5 hdf5 h5py keyutils xz \
      ncurses libnsl readline sqlite openssl libuuid \
      libxml2 libxslt blast r-base swig tzdata wget zlib \
 && conda clean -afy

WORKDIR /app
COPY . /app



# Run project‑specific setup (it will now see the iggypop env)
RUN touch /.dockerenv \
 && chmod +x setup.sh \
 && ./setup.sh \
 && rm /.dockerenv

# Make the entry‑point executable
RUN chmod +x iggypop.py

ENTRYPOINT ["/bin/bash"]
