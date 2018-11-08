#!bin/bash
#
# 
# Usage :
#
# export CI_PROJECT_DIR=<PATH TO SICONOS TUTORIALS REPOSITORY> (or any 'reference path')
# sh install_oce clone_oce
# OR
# sh install_oce
#
# Remark : when call from gitlab ci script, the export CI_PROJECT_DIR step is useless.
#
#
# Result :
# - clone oce (from github) in $CI_PROJECT_DIR/build/
# - configure and build oce in $CI_PROJECT_DIR/build/oce-last
# - install oce in CI_PROJECT_DIR/install/oce
# - clone pythonocc (from github) in $CI_PROJECT_DIR/build/
# - configure and build pythonocc in $CI_PROJECT_DIR/build/pythonocc
# - install pythonocc package in $CI_PROJECT_DIR/install/site-packages
# - set PYTHONPATH properly to allow pythonocc usage (--> test python -c 'import OCC')
#
# When called without 'clone_oce' arg, the first three steps are ignored (install of oce).
# In that case pythonocc is built and installed assuming a former installation of oce (e.g. with package manager)
#

# Get number of procs
if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi
# Check if CI_PROJECT_DIR is set AND not empty
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos-tutorials' repository (absolute) path."}
# Create directory to install libs.
mkdir -p $CI_PROJECT_DIR/install/

# Build dir
mkdir -p build
cd build
export build_path=$PWD

#  -- OCE --
# If 'clone_oce' is set as input arg of the present script,
# download and install the last version of oce.

if [ "$#" -eq 1 ] && [ "$1" = "clone_oce" ]; then 
    # If "clone_oce" as third arg of this script ...
    echo "Clone last version of oce ..."
    git clone https://github.com/tpaviot/oce.git
    mkdir oce-last
    cd oce-last
    # Warning : install in 'user' path, that will be transfered between jobs (artifacts)
    cmake ../oce -DOCE_INSTALL_PREFIX=$CI_PROJECT_DIR/install/oce  -Wno-deprecated -Wno-dev -DCMAKE_BUILD_TYPE=Release
    make -j $nbprocs > /dev/null
    make install > /dev/null
    # Save path to OCEConfig.cmake, required to configure pythonocc
    export OCE_INSTALL=`grep OCEConfig.cmake install_manifest.txt| sed 's/OCEConfig.cmake//g'`
    oce_option="-DOCE_DIR=$OCE_INSTALL"
    cd $build_path
else
    echo "Use installed version of oce (package?)"
    oce_option=""
fi

# -- Python occ --
# Clone last pythonocc version.
# We assume it is complient with the installed oce version.
# Maybe we should clone specific tags for oce and pythonocc? 
git clone https://github.com/tpaviot/pythonocc-core.git
mkdir pythonocc
cd pythonocc
# Requires (in calling script):
# installpath=`python3 -c "import site;print(site.USER_SITE)"`# Unfortunately, this cannot work, artifacts must be
# in CI_PROJECT_DIR ...
export pyocc_installpath=$CI_PROJECT_DIR/install/site-packages
# Mind the OCC at the end of the install path!
cmake ../pythonocc-core -DCMAKE_BUILD_TYPE=Release -Wno-deprecated $oce_option -DPYTHONOCC_INSTALL_DIRECTORY=$pyocc_installpath/OCC
make install -j $nbprocs > /dev/null
cd $build__path
# test ...
export PYTHONPATH=$pyocc_installpath
python3 -c 'import OCC; print(OCC.__file__)'
