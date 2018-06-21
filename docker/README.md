
To build the Docker image, move to the folder with the desired configuration and run Docker build:

    cd incremental_ubuntu_1604\siemens_gadgetron
    docker build --no-cache --force-rm -t siemens_gadgetron .

List docker images
    sudo docker images

List docker containers
    sudo docker ps -a

Save docker image
    sudo docker save -o ~/tmp/siemens_gadgetron.tar [image name]
