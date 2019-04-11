# Pull image from Siconos project registry (with all required deps for Siconos)
FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04
WORKDIR /home
# Clone Siconos
RUN git clone https://github.com/siconos/siconos.git
RUN cd siconos
RUN git rev-parse --short HEAD > $HOME/siconos-commit-number.txt
RUN mkdir -p /home/build; cd /home/build
# Configure (default), build and install (system wide)
RUN cmake /home/siconos -DCMAKE_INSTALL_PREFIX=/home/install-siconos
RUN make -j 4
RUN make install
ENV SICONOS_INSTALL_DIR /home/install-siconos
