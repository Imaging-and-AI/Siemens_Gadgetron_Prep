# Shared arguments
ARG USERNAME="gtuser"
ARG USER_UID=1007
ARG USER_GID=$USER_UID

ARG GADGETRON_BRANCH=TEP_and_OR

ARG INSTALL_HOME=/opt/packages
# ARG INSTALL_HOME=/opt/conda/envs/gadgetron

#ARG BASE_IMAGE=nvidia/cuda:12.1.1-runtime-ubuntu22.04
#ARG BASE_IMAGE=nvidia/cuda:12.6.2-cudnn-devel-ubuntu24.04
#ARG BASE_IMAGE=nvidia/cuda:12.4.1-runtime-ubuntu22.04
#ARG BASE_IMAGE=gadgetronnhlbi/base:nvidia550-ubuntu24.04
# ARG BASE_IMAGE=nvidia/cuda:12.6.3-runtime-ubuntu24.04
# ARG BASE_IMAGE=nvidia/cuda:13.0.2-cudnn-runtime-ubuntu24.04
# ARG BASE_IMAGE=nvidia/cuda:12.8.1-runtime-ubuntu24.04
ARG BASE_IMAGE=ubuntu:resolute

ARG ENV_NAME=gadgetron

#===========================================================================================
# base image
FROM $BASE_IMAGE AS gadgetron_qperf_baseimage

ARG USERNAME
ARG USER_UID
ARG USER_GID
ARG HOME=/home/$USERNAME
ARG INSTALL_HOME
ARG ENV_NAME
ARG OPENRECON_LABEL

RUN apt-get -o "Acquire::https::Verify-Peer=false" update \
#    && DEBIAN_FRONTEND=noninteractive apt-get -o "Acquire::https::Verify-Peer=false" install -y sudo wget git-core rsync curl net-tools libxml2 libsm6 libxext6 libgl1 libglx-mesa0 supervisor jove dos2unix \
&& DEBIAN_FRONTEND=noninteractive apt-get -o "Acquire::https::Verify-Peer=false" install -y sudo wget git-core rsync curl net-tools libsm6 libxext6 libgl1 libglx-mesa0 supervisor jove dos2unix \
&& apt-get -o "Acquire::https::Verify-Peer=false" clean

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    #
    # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

ARG MAMBAFORGE_VERSION=22.9.0-2
ARG CONDA_GID=900

# Based on https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile
RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MAMBAFORGE_VERSION}/Mambaforge-${MAMBAFORGE_VERSION}-Linux-$(uname -m).sh -O /tmp/miniforge.sh \
    && /bin/bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh \
    && /opt/conda/bin/mamba clean --tarballs --index-cache --packages --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete \
    && /opt/conda/bin/mamba clean --force-pkgs-dirs --all --yes  \
    && groupadd -r conda --gid ${CONDA_GID} \
    && usermod -aG conda ${USERNAME} \
    && chown -R $USER_UID:$USER_GID /opt/conda \
    && chmod -R g+w /opt/conda \
    && find /opt -type d | xargs -n 1 chmod g+s

# Copy environment, which will be filtered for later staged
COPY --chown=$USER_UID:$USER_GID ./gadgetron/environment.yml /tmp/build/

# Create mount points for tests
RUN mkdir -p /test && chown ${USER_UID}:${USER_GID} /test && mkdir -p ${INSTALL_HOME}
VOLUME /test

# Add a section to /etc/bash.bashrc that ensures that a section is present at the end of ~/.bashrc.
# We can't just write to .bashrc from here because it will be overwritten if the vscode user has
# opted to use their own dotfiles repo. The dotfiles repo is cloned after the postCreateCommand
# in the devcontainer.json file is executed.
RUN echo "\n\
if ! grep -q \"^source /opt/conda/etc/profile.d/conda.sh\" ${HOME}/.bashrc; then\n\
	echo \"source /opt/conda/etc/profile.d/conda.sh\" >> ${HOME}/.bashrc\n\
	echo \"conda activate $(grep 'name:' /tmp/build/environment.yml | tr -d '\r' | awk '{print $2}')\" >> ${HOME}/.bashrc\n\
fi\n" >> /etc/bash.bashrc

ENV TINI_VERSION=v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

RUN echo 'alias gt_log="tail -F -n2000 /tmp/gadgetron.log"' >> ${HOME}/.bashrc
RUN echo 'alias gt_config="cd /opt/conda/envs/gadgetron-qperf/share/gadgetron/config"' >> ${HOME}/.bashrc
RUN echo 'alias gt_ai_log="tail -F -n2000 /tmp/gadgetron-model.log"' >> ${HOME}/.bashrc

LABEL "com.siemens-healthineers.magneticresonance.openrecon.metadata:1.1.0"=${OPENRECON_LABEL}

ENV GADGETRON_HOME=/opt/conda/envs/${ENV_NAME}
ENV ISMRMRD_HOME=${GADGETRON_HOME}

ENV PATH=$PATH:/opt/conda/condabin:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/gadgetron-qperf/lib/python3.10/site-packages/nvidia/nvjitlink/lib:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV OMP_MAX_ACTIVE_LEVELS=1
ENV MKL_NUM_THREADS=32
ENV MKL_DYNAMIC=true
ENV MKL_CBWR=COMPATIBLE
ENV OMP_WAIT_POLICY=PASSIVE
ENV OPENBLAS_NUM_THREADS=1 \
    LC_ALL=C \
    GTBABYLON_PUBLIC_KEY=/opt/key/gtbabylon-key.public \
    HDF5_USE_FILE_LOCKING=FALSE \
    HDF5_DISABLE_VERSION_CHECK=1 \
    GADGETRON_ISMRMRD_DUMP=OFF

#===========================================================================================
# build conda env

FROM gadgetron_qperf_baseimage AS gadgetron_qperf_base

ARG USERNAME
ARG USER_UID
ARG USER_GID
ARG HOME=/home/$USERNAME
ARG INSTALL_HOME
ARG ENV_NAME
ARG OPENRECON_LABEL

USER ${USER_UID}

RUN mkdir -p ${HOME}/.cache/conda/notices && sudo chown -R ${USER_UID}:$USER_GID ${HOME}/.cache/conda/notices
RUN sudo chown -R $USER_UID:$USER_GID /opt && mkdir -p /opt/code
RUN sudo mkdir -p /opt/integration-test && sudo chown $USER_UID:$USER_GID /opt/integration-test

RUN grep -v "#.*\<cuda\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml && cat /tmp/build/filtered_environment.yml

RUN umask 0002 && /opt/conda/bin/mamba env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/mamba clean -afy && sudo chown -R $USER_UID:$USER_GID /opt/conda && pip3 cache purge
    # . /opt/conda/etc/profile.d/conda.sh && umask 0002 && conda activate ${ENV_NAME} && sh -x && \
    # pip3 uninstall -y torch torchvision && pip3 install --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host download.pytorch.org --pre torch torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cu126 && pip3 cache purge

#===========================================================================================
# build gadgetron

FROM gadgetron_qperf_base AS gadgetron_build

ARG USER_UID
ARG USER_GID
ARG INSTALL_HOME

ARG GADGETRON_BRANCH
ARG ENV_NAME
ARG OPENRECON_LABEL

LABEL "com.siemens-healthineers.magneticresonance.openrecon.metadata:1.1.0"=${OPENRECON_LABEL}

USER ${USER_UID}
WORKDIR /opt
SHELL ["/bin/bash", "-c"]

# RUN . /opt/conda/etc/profile.d/conda.sh && umask 0002 && conda activate ${ENV_NAME} && sh -x && \
#     conda install git

# RUN cd /opt/code && git clone https://github.com/gadgetron/gadgetron --branch ${GADGETRON_BRANCH} --single-branch && cd /opt/code/gadgetron

RUN . /opt/conda/etc/profile.d/conda.sh && umask 0002 && conda activate ${ENV_NAME} && sh -x && \
    cd /opt/code && git clone https://github.com/gadgetron/gadgetron --branch ${GADGETRON_BRANCH} --single-branch && cd /opt/code/gadgetron

# Patch some files to reduce verbosity
COPY --chown=$USER_UID:$USER_GID /patch /opt/code/patch/
RUN cp /opt/code/patch/IsmrmrdDumpGadget.cpp    /opt/code/gadgetron/gadgets/mri_core/IsmrmrdDumpGadget.cpp
RUN cp /opt/code/patch/NoiseAdjustGadget.cpp    /opt/code/gadgetron/gadgets/mri_core/NoiseAdjustGadget.cpp

RUN . /opt/conda/etc/profile.d/conda.sh && umask 0002 && conda activate ${ENV_NAME} && sh -x && \
    cd /opt/code/gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ -DCMAKE_BUILD_TYPE=Release -DUSE_MKL=ON -DUSE_CUDA=OFF -DDISABLE_FORK=OFF -DREQUIRE_SIGNED_CONFIG=OFF -DCMAKE_INSTALL_PREFIX=$INSTALL_HOME && \
    make -j $(nproc) && \
    make install

COPY /Siemens_Gadgetron_Prep /opt/code/Siemens_Gadgetron_Prep/
USER root
RUN chown -R $USER_UID:$USER_GID /opt/code/Siemens_Gadgetron_Prep
USER ${USER_UID}

RUN . /opt/conda/etc/profile.d/conda.sh && umask 0002 && conda activate ${ENV_NAME} && sh -x && \
    cd /opt/code/Siemens_Gadgetron_Prep && \
    mkdir build && \
    cd build && \
    cmake ../ -DCMAKE_BUILD_TYPE=Release -DUSE_MKL=ON -DUSE_CUDA=OFF -DDISABLE_FORK=OFF -DREQUIRE_SIGNED_CONFIG=OFF -DGADGETRON_HOME=$INSTALL_HOME -DCMAKE_INSTALL_PREFIX=$INSTALL_HOME -DENV_NAME=${ENV_NAME} && \
    make -j $(nproc) && \
    make install 
    # git rev-parse HEAD >> /opt/code/Siemens_Gadgetron_Prep_sha1.txt && \
    # /opt/code/gadgetron/docker/manifest --key .io.gadgetron.Siemens_Gadgetron_Prep.sha1 --value `git rev-parse HEAD`

# # Cleanup files not required after compiling
# RUN  rm -r /root/.cache/pip

# # Clean up code
# RUN rm -rf /opt/code/Siemens_Gadgetron_Prep

# ENTRYPOINT [ "/opt/packages/bin/gadgetron" ]

RUN cp /opt/code/gadgetron/docker/entrypoint.sh /opt/
RUN chmod +x /opt/entrypoint.sh
# RUN cp /opt/code/gtprep/docker/conda/start_supervisor /opt/ && chmod +x /opt/start_supervisor
# RUN cp /opt/code/gtprep/docker/conda/supervisord.conf /opt/ && chmod +x /opt/supervisord.conf
# RUN cp -r /opt/code/gadgetron/test/integration /opt/integration-test/
# RUN cp /opt/code/gtprep/test/integration/test_cases.txt /opt/integration-test/test_cases.txt
# RUN cp -r $INSTALL_HOME/* $GADGETRON_HOME/

# RUN rm -rf /opt/code/gadgetron && \
#     rm -rf /opt/code/gtprep && \
#     rm -rf /opt/code/Siemens_Gadgetron && \
#     rm -rf /opt/code/cmr_ml && \
#     rm -rf /opt/code/gt-sim

#CMD ["/opt/conda/bin/conda", "run", "-n", ${ENV_NAME}, "--no-capture-output", "${ENV_NAME}/bin/gadgetron", ">>", "/tmp/gadgetron.log"]
# ENTRYPOINT [ "/tini", "--", "/opt/entrypoint.sh" ]



FROM gadgetron_qperf_baseimage AS gadgetron_gtprep_runtime_service_build

ARG USER_UID
ARG USER_GID
ARG ENV_NAME
ARG OPENRECON_LABEL

ARG USERNAME
ARG HOME=/home/$USERNAME
ARG INSTALL_HOME

USER ${USER_UID}

# LABEL "com.siemens-healthineers.magneticresonance.openrecon.metadata:1.1.0"=${OPENRECON_LABEL}
COPY --from=gadgetron_build --chown=$USER_UID:$USER_GID opt/code/gadgetron/environment.yml /tmp/build/

RUN mkdir -p ${HOME}/.cache/conda/notices && sudo chown -R ${USER_UID}:$USER_GID ${HOME}/.cache/conda/notices
RUN sudo chown -R $USER_UID:$USER_GID /opt && mkdir -p /opt/code
RUN grep -v -vE "#.*\<dev\>|#.*cuda" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml && cat /tmp/build/filtered_environment.yml
RUN umask 0002 && /opt/conda/bin/mamba env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/mamba clean -afy && sudo chown -R $USER_UID:$USER_GID /opt/conda
RUN . /opt/conda/etc/profile.d/conda.sh && umask 0002 && conda activate ${ENV_NAME}&& sh -x && pip3 cache purge
RUN rm -r /tmp/build
# RUN sudo mkdir -p /opt/integration-test && sudo chown $USER_UID:$USER_GID /opt/integration-test
# RUN mkdir -p /opt/conda/envs/gadgetron-gtprep/log/supervisor
# RUN mkdir -p /opt/conda/envs/gadgetron-gtprep/run/supervisor
RUN mkdir -p /tmp/gadgetron_data && chmod 777 /tmp/gadgetron_data
RUN mkdir -p /tmp/gadgetron && chmod 777 /tmp/gadgetron

COPY --from=gadgetron_build --chown=$USER_UID:$USER_GID ${INSTALL_HOME} /opt/conda/envs/${ENV_NAME}
COPY --from=gadgetron_build --chown=$USER_UID:$USER_GID /opt/entrypoint.sh /opt/entrypoint.sh

ENTRYPOINT [ "/tini", "--", "/opt/entrypoint.sh" ]
