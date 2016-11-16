# music-osc
Experimental interface between MUSIC and OSC

Set correct SRV_IP in all source files. Default is localhost.

To build mpi binaries:

  make mpi

To build server binaries:

  make server

or on OS X

  make CXXFLAGS="-std=c++11 -g -O3" server

To make test example:

  make test

Running test example:

On server, launch both udptoosc and osctoudp:

  ./usptoosc & ./osctoudp

Launch test.music:

  mpirun -np 3 music test.music

Output should look like:

udptoosc: 1.004 KEY 3 PRESSED
udptoosc: 1.004 KEY 4 PRESSED
udptoosc: 1.504 KEY 3 RELEASED
udptoosc: 1.504 KEY 4 RELEASED
udptoosc: 2.505 KEY 6 PRESSED
udptoosc: 2.505 KEY 7 PRESSED
udptoosc: 3.005 KEY 6 RELEASED
udptoosc: 3.005 KEY 7 RELEASED
.
.
.
