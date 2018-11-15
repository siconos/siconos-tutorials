FROM base:archlinux
RUN pacman -Sy && pacman --noconfirm -S \
        git \
        boost \
        boost-libs \
        gmp \
        cmake \
        make \
        gcc \
        gcc-fortran \
        swig \
        openblas \
        lapacke \
        cblas \
        python3 \
        python-pytest \
        python-scipy \
        python-h5py \
        python-pip \
        cppunit \
        bullet \
        freetype2 \
        freeglut \
        glu
RUN git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git
ENV CI_PROJECT_DIR /siconos-tutorials
RUN sh $CI_PROJECT_DIR/ci/install_oce clone_oce
RUN cd $CI_PROJECT_DIR/
RUN git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git
RUN mkdir -p build-siconos && cd build-siconos


