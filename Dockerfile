FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"
USER docker
RUN cd /home/docker \
    && wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/docker/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

USER root
# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev && \
    rm -rf /var/lib/apt/lists/*

# Create Conda environments
USER docker
RUN conda install -c bioconda -c conda-forge nanoplot=1.40.0
RUN conda install -c bioconda -c conda-forge nanofilt=2.8.0
RUN conda install -c bioconda -c conda-forge flye
RUN conda config --set solver classic
RUN conda install -c bioconda -c conda-forge minimap2
RUN conda install -c bioconda -c conda-forge racon=1.5.0
RUN conda install -c bioconda -c conda-forge checkm-genome
RUN conda install -c bioconda -c conda-forge bwa
RUN conda install -c bioconda -c conda-forge samtools
RUN conda install -c bioconda -c conda-forge pilon

RUN conda create -n medaka -c conda-forge -c bioconda medaka

USER root
RUN apt-get update && apt-get install -y pkg-config libfreetype6-dev libpng-dev python3-matplotlib
USER docker
RUN cd /home/docker/ && git clone https://github.com/ablab/quast && \ 
    cd quast && ./setup.py install_full

#RUN conda install -c bioconda -c conda-forge quast=5.2.0

# Install CheckM database
RUN echo "export CHECKM_DATA_PATH=/home/docker/checkm" >> /home/docker/.bashrc
ENV CHECKM_DATA_PATH=/home/docker/checkm
RUN mkdir /home/docker/checkm && cd /home/docker/checkm && \
    wget -q https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar -xvzf checkm_data_2015_01_16.tar.gz && rm checkm_data_2015_01_16.tar.gz
RUN mkdir /home/docker/code
COPY ONTAssenbler.sh /home/docker/code/ONTAssenbler.sh
#RUN chmod +x /home/docker/code/ONTAssenbler.sh

# Set working directory
WORKDIR /Data
