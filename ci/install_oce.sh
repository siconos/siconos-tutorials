#!bin/bash
#

# Get number of procs
if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi

mkdir -p build
cd build
export build_path=$PWD

#  -- OCE --

if [ "$#" -eq 1 ] && [ "$1" = "clone_oce" ]; then 
    # If "clone_oce" as third arg of this script ...
    echo "Clone last version of oce ..."
    git clone https://github.com/tpaviot/oce.git
    mkdir oce-last
    cd oce-last
    cmake ../oce -DOCE_INSTALL_PREFIX=$HOME/install_oce  -Wno-deprecated -DCMAKE_BUILD_TYPE=Release
    make -j $nbprocs > /dev/null
    make install > /dev/null
    cd $build_path
    export OCE_INSTALL=$HOME/install_oce
else
    echo "Use installed version of oce (package?)"
fi

# -- Python occ --
git clone https://github.com/tpaviot/pythonocc-core.git
mkdir pythonocc
cd pythonocc
# Requires (in calling script):
# export installpath=`python3 -c "import site;print(site.USER_SITE)"`
cmake ../pythonocc-core -DCMAKE_BUILD_TYPE=Release -Wno-deprecated -DOCE_DIR=$OCE_INSTALL -DPYTHONOCC_INSTALL_DIRECTORY=$installpath
make install -j $nbprocs > /dev/null
cd $build__path
# test ...
python3 -c 'import OCC; print(OCC.__file__)'
