FROM condaforge/mambaforge

ENV INSTALL_DIR="/opt"

RUN apt-get update; export DEBIAN_FRONTEND=noninteractive; apt-get -y --no-install-recommends install \
   wget \
   unzip \
   build-essential \
   gcc-multilib \
   multiqc \
   netcat-openbsd && \
   export DEBIAN_FRONTEND=noninteractive; apt-get -y install openssh-client openssh-sftp-server openssh-server

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH
# Install the conda environments
COPY * /
RUN sh installEnvs.sh
SHELL ["/bin/bash", "--login", "-c"]
RUN source /opt/conda/etc/profile.d/conda.sh
USER 1000
#RUN source activate
