#!bin/bash
#
export build_path=$PWD

#  -- OCE --

if [ "$#" -eq 1 ] && [ "$1" = "clone_oce" ]; then 
    # If "clone_oce" as third arg of this script ...
    echo "Clone last version of oce ..."
    git clone https://github.com/tpaviot/oce.git
    mkdir oce-last
    cd oce-last
    cmake ../../oce -DOCE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release
    make -j $nbprocs
    make install
    cd $build_path
else
    echo "Use installed version of oce (package?)"
fi

# -- Python occ --
git clone https://github.com/tpaviot/pythonocc-core.git
mkdir pythonocc
cd pythonocc
cmake ../pythonocc-core -DCMAKE_BUILD_TYPE=Release
make install -j $nbprocs
cd $build__path
# test ...
python3 -c 'import OCC'
