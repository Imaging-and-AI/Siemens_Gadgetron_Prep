FROM gadgetron/ubuntu_2004

RUN rm -rf /opt/code/gadgetron

RUN cd /opt/code && \
    git clone https://github.com/gadgetron/gadgetron.git && \
    cd gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` && \
    cp /opt/code/gadgetron/docker/start_supervisor /opt/ && \
    cp /opt/code/gadgetron/docker/supervisord.conf /opt/

COPY /Siemens_Gadgetron_Prep-T1_4.1 /opt/code/Siemens_Gadgetron_Prep/

RUN cd /opt/code/Siemens_Gadgetron_Prep && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j 8 && \
    make install 
    # git rev-parse HEAD >> /opt/code/Siemens_Gadgetron_Prep_sha1.txt && \
    # /opt/code/gadgetron/docker/manifest --key .io.gadgetron.Siemens_Gadgetron_Prep.sha1 --value `git rev-parse HEAD`

RUN pip3 install antspyx dill multiprocessing_on_dill

# Clean up code
RUN rm -rf /opt/code/Siemens_Gadgetron_Prep

CMD ["/usr/local/bin/gadgetron"]
