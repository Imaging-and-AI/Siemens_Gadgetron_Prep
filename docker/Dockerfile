FROM ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_2004:latest

RUN  mkdir -p /opt/code

# --- Re-compile Gadgetron if needed ---
# RUN rm -rf /opt/code/gadgetron

# RUN cd /opt/code && \
#     git clone https://github.com/gadgetron/gadgetron.git

# # Copy in locally modified files
# COPY /Server.cpp            /opt/code/gadgetron/apps/gadgetron/
# COPY /External.cpp          /opt/code/gadgetron/apps/gadgetron/connection/stream/common/
# COPY /stream/External.cpp   /opt/code/gadgetron/apps/gadgetron/connection/stream/
# COPY /start-gadgetron.sh.in /opt/code/gadgetron/chroot/

# # Check out specific branch
# RUN cd /opt/code/gadgetron && git fetch origin pull/880/head:IPv6fix && git checkout IPv6fix

# RUN cd /opt/code/gadgetron && \
#     mkdir build && \
#     cd build && \
#     cmake ../ && \
#     make -j $(nproc) && \
#     make install && \
#     /opt/code/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` && \
#     cp /opt/code/gadgetron/docker/start_supervisor /opt/ && \
#     cp /opt/code/gadgetron/docker/supervisord.conf /opt/


# siemens_to_ismrmrd
RUN rm -rf /opt/code/siemens_to_ismrmrd
RUN cd /opt/code && \
    git clone https://github.com/kspaceKelvin/siemens_to_ismrmrd.git && \
    cd /opt/code/siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

# ANTs dependencies for mSASHA
RUN pip3 install --no-cache-dir antspyx dill multiprocessing_on_dill

COPY /Siemens_Gadgetron_Prep /opt/code/Siemens_Gadgetron_Prep/

RUN cd /opt/code/Siemens_Gadgetron_Prep && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j 8 && \
    make install 
    # git rev-parse HEAD >> /opt/code/Siemens_Gadgetron_Prep_sha1.txt && \
    # /opt/code/gadgetron/docker/manifest --key .io.gadgetron.Siemens_Gadgetron_Prep.sha1 --value `git rev-parse HEAD`

# Cleanup files not required after compiling
RUN  rm -r /root/.cache/pip

# Clean up code
RUN rm -rf /opt/code/Siemens_Gadgetron_Prep

CMD ["/usr/local/bin/gadgetron"]
