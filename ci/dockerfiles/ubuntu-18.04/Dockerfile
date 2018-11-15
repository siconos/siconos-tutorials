FROM ubuntu:18.04
RUN  apt update && apt install -y -qq \
        cmake \
        git-core \
        make \
        libboost-dev \
        libgmp-dev \
        swig \
        gcc \
        gfortran \
        g++ \
        liblapack-dev \
        libatlas-base-dev \
        lp-solve \
        liblpsolve55-dev \
        libpython3-dev \
        doxygen \
        libcppunit-dev \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev \
        python3-pip\
        python3-scipy \
        python3-pytest \
        valgrind \
        python3-h5py \
        liboce-foundation-dev \
        liboce-modeling-dev \
        liboce-ocaf-dev \
        liboce-visualization-dev
RUN mkdir -p /builds
ENV CI_PROJECT_DIR /builds
RUN cd $CI_PROJECT_DIR && git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git
RUN cd $CI_PROJECT_DIR && git clone https://github.com/siconos/siconos.git
RUN cd $CI_PROJECT_DIR && mkdir build-siconos
RUN cd $CI_PROJECT_DIR/build-siconos && cmake -DUSER_OPTIONS_FILE=${CI_PROJECT_DIR}/siconos-tutorials/ci/siconos_with_mechanisms.cmake -DCMAKE_CXX_STANDARD=11 -DSICONOS_USE_BOOST_FOR_CXX11=OFF $CI_PROJECT_DIR/siconos \
        && make
RUN cd $CI_PROJECT_DIR/build-siconos && make install


