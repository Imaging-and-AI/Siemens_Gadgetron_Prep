#FROM gadgetron/ubuntu_1804_base_basic
#FROM kspacekelvin/gadgetron

FROM kspacekelvin/siemens_gtprep_intermediate_1804:latest

# Fix missing copy
RUN cp -f /opt/code/gadgetron/core/variant.hpp                          /usr/local/include/variant.hpp
RUN cp -f /opt/code/gadgetron/core/MessageID.h                          /usr/local/include/MessageID.h
RUN cp -f /opt/code/gadgetron/core/io/from_string.hpp                   /usr/local/include/from_string.hpp
RUN cp -f /opt/code/gadgetron/gadgets/mri_core/ImageArraySendMixin.h    /usr/local/include/ImageArraySendMixin.h
RUN cp -f /opt/code/gadgetron/gadgets/mri_core/ImageArraySendMixin.hpp  /usr/local/include/ImageArraySendMixin.hpp

# New version to support changing port number
COPY /start-gadgetron.sh /usr/local/share/gadgetron/chroot/

#Assume that the context is the top level of the git repo
COPY /Siemens_Gadgetron_Prep /opt/code/Siemens_Gadgetron_Prep/

# Remove build folder, if exists
RUN rm -rf /opt/code/Siemens_Gadgetron_Prep/build

#gtprep compilation
RUN cd /opt/code/Siemens_Gadgetron_Prep && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j 8 && \
    git rev-parse HEAD >> /opt/code/Siemens_Gadgetron_Prep_sha1.txt && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.Siemens_Gadgetron_Prep.sha1 --value `git rev-parse HEAD` && \
    make install 

#Copy test case file
RUN cp -f /opt/code/Siemens_Gadgetron_Prep/test/integration/test_cases.txt /opt/code/test_cases.txt

# #Clean up code
# RUN rm -rf /opt/code/gadgetron
# RUN rm -rf /opt/code/ismrmrd
# RUN rm -rf /opt/code/siemens_to_ismrmrd
# RUN rm -rf /opt/code/ismrmrd-python
# RUN rm -rf /opt/code/ismrmrd-python-tools
# RUN rm -rf /opt/code/Siemens_Gadgetron_Prep

# #Clean some MKL stuff
# RUN rm -rf /opt/intel/mkl/include

# # Clean up packages.
# RUN  apt-get clean && \
#    rm -rf /var/lib/apt/lists/*

# # Remove Python files to reduce size
# RUN pip  freeze | xargs -n 1 pip  uninstall -y; exit 0
# RUN pip3 freeze | xargs -n 1 pip3 uninstall -y; exit 0
# RUN rm -rf /root/.cache/pip

CMD ["/usr/local/bin/gadgetron"]
