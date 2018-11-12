Siemens_Gadgetron Docker Image Generation
------------------------------

Since the Siemens_Gadgetron code cannot be cloned inside the container (it would need git access), we will assume that the docker build context is at the root of the Siemens_Gadgetron repo. 
The Siemens_Gadgetron repo will be added into the docker image and built after gadgetron.

```
cd ~/mprogs
git clone git@github.com:NHLBI-MR/Siemens_Gadgetron_Prep.git
docker build --network=host -t siemens_gadgetron_prep/gadgetron_siemens_ubuntu1604_cuda90 -f Siemens_Gadgetron_Prep/docker/incremental_ubuntu_1604/siemens_gadgetron_cuda90/Dockerfile .
```

And done. 

This Docker recipe is an incremental recipe from the main Gadgetron repository. 
