# Base image
FROM google/deepvariant:1.5.0

# Default work directory
WORKDIR /usr/src/app
# Avoid interactive prompts for tzdata
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies and AWS CLI
RUN apt-get update && apt-get install -y \
    software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && apt-get install -y \
    bash \
    curl \
    tzdata \
    wget \
    unzip \
    bzip2 \
    gcc \
    make \
    samtools \
    bwa \
    awscli \
    wget \
    bcftools \
    nano \
    less \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    pkg-config \
    zlib1g-dev \
    perl \
    cpanminus \
    libmysqlclient-dev \
    libxml2 \
    libexpat1-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    zlib1g-dev  \
    tabix \
    python3 \
    python3.8 \
    python3.8-distutils \
    python3.10 \
    python3.10-distutils \
    python3-pip && \
    ln -fs /usr/share/zoneinfo/Etc/UTC /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda for Conda environment management
RUN mkdir -p /opt/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda3/miniconda.sh && \
    bash /opt/miniconda3/miniconda.sh -b -u -p /opt/miniconda3 && \
    rm /opt/miniconda3/miniconda.sh
ENV PATH="/opt/miniconda3/bin:${PATH}"
# Verify conda installation
RUN conda --version

RUN conda config --set solver classic
# Disable plugins during Conda update to bypass errors
RUN CONDA_NO_PLUGINS=true conda update -n base -c defaults conda -y

# Install OpenJDK 17
RUN apt-get update && apt-get install -y openjdk-17-jre-headless && apt-get clean
# Set JAVA_HOME environment variable
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH="$JAVA_HOME/bin:${PATH}"

# Install Python 3.8 and its dependencies
RUN python3.8 -m pip install --upgrade pip && \
    python3.8 -m pip install tensorflow==2.11.0 absl-py protobuf

# Create a wrapper script for DeepVariant
# Fix the run_deepvariant script to use python3.8
RUN echo '#!/bin/bash\n\
exec python3.8 -u /opt/deepvariant/bin/run_deepvariant.py "$@"' > /opt/deepvariant/bin/run_deepvariant && \
    chmod +x /opt/deepvariant/bin/run_deepvariant

# Install FastQC
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -O /tmp/fastqc.zip && \
    unzip /tmp/fastqc.zip -d /opt && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /tmp/fastqc.zip

# Install Fastp (specific version)
RUN wget -q http://opengene.org/fastp/fastp.0.23.4 -O /tmp/fastp && \
    chmod +x /tmp/fastp && \
    mv /tmp/fastp /usr/local/bin/

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ \
    && apt-get update && apt-get install -y gpg graphviz

# Install Docker CLI inside the container
RUN apt-get update && apt-get install -y \
    docker.io \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.27.1/picard.jar -O /usr/local/bin/picard.jar \
    && chmod +x /usr/local/bin/picard.jar

# Install perl module
RUN apt-get update && \
    apt-get -y install gcc perl cpanminus libmysqlclient-dev libxml2 libexpat1-dev libbz2-dev liblzma-dev libxml2-dev zlib1g-dev && \
    apt-get clean

# Install the required Perl modules before installing BioPerl
RUN cpanm --force Net::HTTP && \
    cpanm --force XML::Parser::PerlSAX && \
    cpanm XML::Twig && \
    cpanm LWP::UserAgent && \
    cpanm XML::DOM && \
    cpanm --force XML::LibXML && \
    cpanm --force XML::LibXML::Reader && \
    cpanm App::cpanminus Config::Tiny && \
    cpanm Module::Build && \
    cpanm Test::Warnings && \
    cpanm DBI && \
    cpanm DBD::mysql && \
    cpanm --force BioPerl && \
    cpanm --force Bio::Perl && \
    apt-get clean

# Clone and install Ensembl API
RUN rm -rf /usr/src/app/ensembl-vep && \
    git clone https://github.com/ensembl/ensembl-vep /usr/src/app/ensembl-vep && \
    cd /usr/src/app/ensembl-vep && \
    cpanm --installdeps . && \
    cpanm Module::Build && \
    perl INSTALL.pl && \
    perl INSTALL.pl -a p --PLUGINS all --PLUGINSDIR ./Plugins/

# REFERENCE FILES COPY FROM s3


# Install MANTA
COPY p2manta.yml /tmp/
RUN conda env create -f /tmp/p2manta.yml

## install Mutserve
RUN mkdir mutserveTool && \
    cd mutserveTool && \
    curl -sL mutserve.vercel.app  | bash

# Download database zip files , exomiser zip file, application.properties
RUN apt-get install -y p7zip-full
RUN mkdir -p /usr/src/app/exomiser/db 
COPY exomiser-cli-14.0.0-distribution.zip /tmp/exomiser-cli-14.0.0-distribution.zip && \
    unzip /tmp/exomiser-cli-14.0.0-distribution.zip -d /usr/src/app/exomiser && \
    rm /tmp/exomiser-cli-14.0.0-distribution.zip
### exomiser database hg37
COPY application.properties /usr/src/app/exomiser/exomiser-cli-14.0.0

# copy exomiser template
RUN mkdir -p /usr/src/app/exomiser/db && \
    cp exomiser-template.yml /usr/src/app/exomiser/db

# install vus prize
COPY vusprize1.yml /tmp/
RUN conda env create -f /tmp/vusprize1.yml
RUN git clone https://github.com/danielhmahecha/VusPrize.git

# Add sitecustomize.py to modify sys.path
# temporary comment
#RUN echo 'import sys\nsys.path = [path for path in sys.path if "/usr/local/lib" not in path]\nprint("Modified sys.path:", sys.path)' \
#    > /opt/miniconda3/envs/vusprize_env/lib/python3.7/site-packages/sitecustomize.py

# Create conda env for vus prize 
# temporary comment
RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && \
    mv bin/micromamba /usr/local/bin/micromamba && \
    rm -rf bin
ENV MAMBA_ROOT_PREFIX=/usr/local/mamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

# Copy pipeline files (optional, if needed during the build process)
#COPY main_icgeb.nf /usr/src/app/main_icgeb.nf
#COPY modules /usr/src/app/modules
#COPY nextflow_updated.config /usr/src/app/nextflow.config
RUN pip install --no-cache-dir numpy pandas openpyxl
# Add a health check
HEALTHCHECK --interval=30s --timeout=10s \
    CMD nextflow -v || exit 1
RUN apt-get update && apt-get install -y jq

# Install R
RUN sudo apt update && \
    sudo apt install software-properties-common dirmngr -y && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    sudo apt update && \
    sudo apt install r-base r-base-dev -y

# Install required R packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(version = '3.21')" && \
    R -e "BiocManager::install(c('TxDb.Hsapiens.UCSC.hg38.knownGene', 'optparse'))"

RUN sudo apt install -y libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev && \
    R -e "BiocManager::install('devtools')"
RUN R -e "devtools::install_github('vplagnol/ExomeDepth')"

RUN cpanm --force Parse::CSV

# Install AnnotSV
RUN git clone https://github.com/lgmgeo/AnnotSV.git /usr/src/app/AnnotSV && \
    cd /usr/src/app/AnnotSV && \
    make PREFIX=. install && \
    make PREFIX=. install-human-annotation
RUN sudo apt install -y bedtools

RUN  wget https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.4.1/glnexus_cli -O /usr/local/bin/glnexus_cli \
    && chmod +x /usr/local/bin/glnexus_cli
# Install BIAS Nirvana
RUN wget https://packages.microsoft.com/config/ubuntu/22.04/packages-microsoft-prod.deb -O packages-microsoft-prod.deb && \
    sudo dpkg -i packages-microsoft-prod.deb && \
    rm packages-microsoft-prod.deb && \
    sudo apt update && \
    sudo apt install -y dotnet-sdk-6.0 && \
    wget https://github.com/Illumina/Nirvana/releases/download/v3.18.1/Nirvana-3.18.1-net6.0.zip -O Nirvana-3.18.1-net6.0.zip && \
    unzip Nirvana-3.18.1-net6.0.zip && \
    rm Nirvana-3.18.1-net6.0.zip && \
    git clone https://github.com/bitscopic/BIAS-2015.git /usr/src/app/BIAS-2015 && \
    conda install -c bioconda sambamba


    # Set entrypoint
#ENTRYPOINT ["/usr/local/bin/fetch_references.sh"]

 # Default command
CMD ["all", "--cores", "5"]
#CMD ["/usr/local/bin/fetch_references.sh", "&&", "tail", "-f", "/dev/null"]
#CMD ["bash", "-c", "tail -f /dev/null"]
# CMD bash -c "\
# #gpg --batch --yes --passphrase '$GPG_PASSPHRASE' -o main_securing3rd.nf -d main_securing3rd.nf.gpg && \
# nextflow -log /usr/src/app/output/pipeline.log run main_icgeb.nf -c nextflow.config -params-file config.json -with-report /usr/src/app/output/pipeline_report -with-dag /usr/src/app/output/pipeline_DAG.png -with-trace /usr/src/app/output/pipeline_trace.txt -with-timeline /usr/src/app/output/pipeline_timeline.html && \
# chmod -R 777 /usr/src/app/output"
RUN conda init bash
