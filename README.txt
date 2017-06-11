++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++                                                      ++
++       WSA: Weakly-hard Schedulability Analyzer       ++
++                                                      ++
++            Author: xxxxxxxxxxxx                      ++
++          xxxxxxxxxxxxxxxxxxxxxxxx                    ++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



bin/            the binary distribution for Linux
src/            source codes
examples/       an example task system
README


Usage
======
Usage:   wsa [-option] [argument]
option:  -h  show help information
         -n number of tasks
         -m  the m parameter
         -k  the K parameter
         -f  the file containing input tasks in the form of (C, D, T)

Example
========
cd wsa/
./bin/wsa -n 6 -m 1 -k 5 -f examples/taskset1/system.dat 



To build
=========
cd wsa/src
make clean
make all
