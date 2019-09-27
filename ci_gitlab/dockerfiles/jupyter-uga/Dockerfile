FROM jupyter/r-notebook
USER root
RUN  apt update && apt install -y -qq \
        cmake \
        libboost-dev \
        libgmp-dev \
        swig \
        liblapack-dev \
        libatlas-base-dev \
        lp-solve \
        liblpsolve55-dev \
        libpython3-dev \
        doxygen \
        libcppunit-dev \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev
RUN conda install --yes scipy \
        h5py \
        pytest
USER $NB_USER
RUN git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git
ENV CI_PROJECT_DIR $HOME/siconos-tutorials
RUN git clone https://github.com/siconos/siconos.git
RUN mkdir $HOME/build-siconos
RUN cd $HOME/build-siconos && cmake -DUSER_OPTIONS_FILE=${CI_PROJECT_DIR}/binder/siconos_confs/siconos_notebook.cmake -DCMAKE_CXX_STANDARD=11 -DSICONOS_USE_BOOST_FOR_CXX11=OFF $HOME/siconos \
        && make
USER root
RUN cd /home/$NB_USER/build-siconos && make install
USER $NB_USER

