# This directory contains the files used to configure a binder repository for siconos.

## About 

Only the Dockerfile is used, the other files are kept for legacy.

The binder image starts from the docker image 'siconoslab-release-4.4' saved in the registries of the project.

It has to be build only once (the result is kept on binder hub) but if the notebooks fails to start, try to re-built it as explained below.

## How to Build the binder interactive siconos notebooks 

Steps to reproduce to create an interactive siconos-tutorials session on binder

* connect to https://mybinder.org/
* 'GitHub repository name or URL' :
    * choose Git repository in the drop-down menu
    * set its content with "https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git"
* Let default value for Git ref
* URL to open (optional) : siconos-tutorials/siconos-notebooks/index.ipynb (to start the session on this file)
* BEFORE launching, copy the given URL (in the field "copy the URL below ...) and update the content of the binder badge in the [readme file](../readme.md) of this project
* Launch and wait ...

