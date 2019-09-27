# Pull image from Siconos project registry (with all required deps for Siconos)
FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/archlinux-oce
WORKDIR /home
# Clone Siconos
RUN git clone https://github.com/siconos/siconos.git
# Warning : from this point we use conf from siconos repos! 
ENV CI_PROJECT_DIR /home/siconos
ENV ctest_build_model Experimental
ENV IMAGE_NAME archlinux-oce
ENV PYTHONPATH /home/install/site-packages
ENV cdash_submit 1
ENV OCE_ROOT /home/install/oce
# RUN sh /home/siconos/ci_gitlab/install_siconos.sh /home/siconos/ci_gitlab/siconos_confs/siconos_with_mechanisms.cmake
# RUN cd /home/siconos/build ; make -j 4 install # install in /home/siconos/install-siconos
ENV SICONOS_INSTALL_DIR /home/siconos/install-siconos
# ENV CI_PROJECT_DIR /home/


