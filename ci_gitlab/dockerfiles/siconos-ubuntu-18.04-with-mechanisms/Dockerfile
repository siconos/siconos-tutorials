# Pull image from Siconos project registry
# (with all required deps for Siconos, including bullet, oce and pythonocc.)
FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04-oce
WORKDIR /home
# Clone Siconos
RUN git clone https://github.com/siconos/siconos.git
# Copy Siconos config into the Docker container
COPY siconos_with_mechanisms.cmake /home/siconos_with_mechanisms.cmake
RUN cd siconos ; git rev-parse --short HEAD > /home/siconos-commit-number.txt
RUN mkdir -p /home/build
# Configure (default), build and install (system wide)
RUN cd /home/build ; cmake /home/siconos -DUSER_OPTIONS_FILE=/home/siconos_with_mechanisms.cmake -DCMAKE_INSTALL_PREFIX=/home/install-siconos
RUN cd /home/build ; make -j 4 ; make install
ENV SICONOS_INSTALL_DIR /home/install-siconos
