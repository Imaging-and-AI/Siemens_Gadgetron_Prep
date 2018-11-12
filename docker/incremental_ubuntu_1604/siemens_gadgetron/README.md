Siemens_Gadgetron Docker Image Generation
------------------------------

Since the Siemens_Gadgetron code cannot be cloned inside the container (it would need git access), we will assume that the docker build context is at the root of the Siemens_Gadgetron repo. 
The Siemens_Gadgetron repo will be added into the docker image and built after gadgetron.

```
git clone git@github.com:NHLBI-MR/Siemens_Gadgetron_Prep.git
cd Siemens_Gadgetron_Prep
docker build -t Siemens_Gadgetron_Prep/gadgetron_siemens_ubuntu1604 -f docker/incremental_ubuntu_1604/siemens_gadgetron/Dockerfile .
```

And done. 

This Docker recipe is an incremental recipe from the main Gadgetron repository. 
