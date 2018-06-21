FROM gadgetron/ubuntu1604_cuda80_cudnn6_base

ARG GADGETRON_URL=https://github.com/gadgetron/gadgetron
ARG GADGETRON_BRANCH=master

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    net-tools cpio

RUN mkdir /opt/code

#ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j 16 && \
    make install

#Fix compiler error with CUDA for now
RUN sed -i "s;#error -- unsupported GNU version! gcc versions later than 4.9 are not supported!;//#error -- unsupported GNU version! gcc versions later than 4.9 are not supported!;g" /usr/local/cuda/include/host_config.h

#GADGETRON
RUN cd /opt/code && \
    git clone ${GADGETRON_URL} --branch ${GADGETRON_BRANCH} --single-branch && \
    cd gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j 16 && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` && \
    cp ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml.example ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml && \
    cp /opt/code/gadgetron/docker/start_supervisor /opt/ && \
    cp /opt/code/gadgetron/docker/supervisord.conf /opt/

#HASH for ISMRMRD
RUN cd /opt/code/ismrmrd && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrd.sha1 --value `git rev-parse HEAD` 

#ISMRMRD PYTHON API
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python.git &&  \
    cd ismrmrd-python && \
    python setup.py install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrdpython.sha1 --value `git rev-parse HEAD` 

#ISMRMRD PYTHON TOOLS
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python-tools.git &&  \
    cd ismrmrd-python-tools && \
    python setup.py install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrdpythontools.sha1 --value `git rev-parse HEAD` 

#SIEMENS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j 16 && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.siemens_to_ismrmrd.sha1 --value `git rev-parse HEAD` 

#Assume that the context is the top level of the git repo
COPY /Siemens_Gadgetron-source /opt/code/Siemens_Gadgetron/

#Siemens_Gadgetron compilation
RUN cd /opt/code/Siemens_Gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j 16 && \
    git rev-parse HEAD >> /opt/code/Siemens_Gadgetron_sha1.txt && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.Siemens_Gadgetron.sha1 --value `git rev-parse HEAD` && \
    make install 

#Copy test case file
RUN cp -f /opt/code/Siemens_Gadgetron/test/integration/test_cases.txt /opt/code/test_cases.txt
    
#Clean up code
RUN rm -rf /opt/code/Siemens_Gadgetron

#Clean some MKL stuff
RUN rm -rf /opt/intel/mkl/include

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*

CMD ["/opt/start_supervisor"]
