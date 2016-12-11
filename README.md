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

  ./udptoosc & ./osctoudp &

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


To compile on Milner, do:

  module swap PrgEnv-cray PrgEnv-gnu
  module add music

  make MPICXX=CC CXX=CC CXXFLAGS="-O3 -std=c++11 -DHAVE_CLOCK_NANOSLEEP"

How to install liblo to $HOME/local, which is where the Makefile expects it

  download liblo-0.28.tar from http://liblo.sourceforge.net
  tar xf liblo-0.28.tar
  cd  liblo-0.28
  ./configure --prefix=$HOME/local --enable-static=yes
  make install


How to run on osc test milner

  on wand: sudo redirect add <IP> 9930 57120
  music -e testosc.music
  chmod +x sim.sh oscin.sh oscout.sh
  cp sim.sh oscin.sh oscout.sh oscin oscout simproxy /cfs/milner/scratch/...
  salloc -N 3
  aprun -n 1 ./sim.sh : -n 1 ./oscin.sh : -n 1 ./oscout.sh
  exit

How to run minimalbrain2 on milner, using filetobrain and braintofile

  The control of the parameters and mode of minimalbrain2 is via the
  parameter file params_3.par which is read during startup. Most
  importantly the mode of the program is set by ctrlmode, currently to
  { 0, 1, 2, 10, 100 }. 10 is the most promising mode at the moment,
  and should allow some kind of interactive runs.

  The version of minimaltest is minimaltest4.music. Use salloc and
  aprun e.g. like:

  salloc -N 6 -t 1:00:00
  aprun -n 1 ./source.sh : -n 1 ./sink.sh : -n 152 ./brain.sh

  Logging can be turned on and off from params_3.par. Note that
  logging slows down a bit.
  
  There are now minitest01.music and minitest02.music that
  demonstrates pattern replay as well as completion ok. NOTE 1 that
  the number of nodes needed is now 6. Also you need to checkout the
  pbcpnn branch ala_changes for now
