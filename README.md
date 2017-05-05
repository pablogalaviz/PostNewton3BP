# olliptic_g4

olliptic_g4 evolves Elliptic-Hyperbolic PDF system   

## Synopsis

Olliptic's current main application is evolutiuon of Schroedinger-Poisson system


## dependencies

## ubuntu 16.04 LTS

Install:

`sudo apt-get install  openmpi-bin libhdf5-openmpi-10-dbg libhdf5-openmpi-10 libhdf5-openmpi-dev libopenmpi-dev libopenmpi1.10 hdf5-tools h5utils libboost-all-dev`

## Mac os x with fink

Install:

`fink install openmpi openmpi-shlibs`


Download boost >= 1.63 and install it at /usr/local/:

`tar -zxvf boost_1_63_0.tar.gz`

`cd boost_1_63_0`

`./b2 install`

Download hdf5 >= 1.8.18 and compile it with openMPI:

`tar -zxvf hdf5-1.8.18.tar.gz`

`cd hdf5-1.8.18`

`CC=/sw/bin/mpicc ./configure --enable-parallel --prefix=/usr/local/`

`make`

`make check`

`sudo make install`


## Compilation

run:

`make`

run 

`make help` to see more options 

## Examples

try

`./Olliptic -h` 
