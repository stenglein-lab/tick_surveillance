## Create a singularity image that includes R libraries needed for this pipeline

This directory contains a singularity definition file to create a singularity image containing the R libraries needed by scripts in this pipeline.  It is based on the [rocker/r-ver docker image](https://hub.docker.com/r/rocker/r-ver) and adds on tidyverse packages and a couple additional packages.  

#### To modify and create the singularity image

Modify the `r_tools.def` file as appropriate

```
# create the singularity image
sudo singularity build r_tools.sif r_tools.def
```

#### To push this image to the singularity library 

I first had to create a singularity library account and login as described [here](https://sylabs.io/guides/latest/user-guide/cloud_library.html?highlight=push#overview)

```
# modify version number as appropriate
singularity push -D "A singularity image containing tidyverse and a few other R packages" r_tools.sif library://stenglein-lab/r_tools/r_tools:1.0.0
```

