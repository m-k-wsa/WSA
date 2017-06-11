++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++                                                      ++
++       WSA: Weakly-hard Schedulability Analyzer       ++
++                                                      ++
++            Author: xxxxxxxxxxxx                      ++
++          xxxxxxxxxxxxxxxxxxxxxxxx                    ++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The weakly hard scheduling analyzes the number of temporal
constraints violations given a time window or a sequence
of task activations. This is also called m-K model, since
a typical formulation consists in checking that no more 
than m deadlines are missed over a set of K activations.
The rationale behind this model is that each deadline miss
will bring the system closer to a faulty condition or
state, and each completion in time will move the system
back towards a safe working state.

In this work, we provide a generalized framework for
offset-free systems scheduled on uniprocessors. Our
formulation is based on a MILP formulation that can be
easily reused to check the system for a large set of m
and K values, allowing the designer to explore a wide
design range. The developed MILP model serves as an
over-approximate analysis, that is, if our weakly hard 
analysis confirms the m-K property, it is guaranteed that
there will be no more than m deadline misses out of any K
successive job activations (of the target task), however,
if the m-K property is not confirmed, we cannot conclude
the opposite.



bin/            the binary distribution for Linux
src/            source codes
examples/       an example task system
README


Usage
======
Usage:   wsa [-option] [argument]
option:  -h  show help information
         -n  number of tasks
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
