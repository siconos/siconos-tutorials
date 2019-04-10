FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04
WORKDIR /home
RUN git clone https://github.com/siconos/siconos.git
RUN mkdir build; cd build
RUN cmake /home/siconos
RUN make -j 4
RUN make install
