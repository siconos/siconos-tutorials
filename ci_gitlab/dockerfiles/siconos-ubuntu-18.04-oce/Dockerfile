FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04-oce
WORKDIR /home
COPY siconos-tutorials-options.cmake .
RUN git clone https://github.com/siconos/siconos.git
RUN mkdir build; cd build
RUN cmake /home/siconos -DUSER_OPTIONS_FILE=/home/siconos-tutorials-options.cmake
RUN make -j 4
RUN make install
