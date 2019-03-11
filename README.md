# DFTBopt: a package to optimize repulsive potentials and electronic parameters for the Density Functional Tight Binding (DFTB) method.

The package includes two main programs: 

(1) repopt for optimization of only repulsive potentials 

(2) erepopt for optimization of all DFTB parameters

erepopt uses repopt for the repulsive potentials fitting. Currently, repopt
is stable released version 1.0. On the other hand, erepopt is still in the beta-development version.
In addition, the package also provide some tools to analyze or evaluate new DFTB parameters.

Requirements
============

In order to install DFTBopt, you need the following software components:

* GNU make

* g++ or icpc compiler

* galib247 and eigen3 library


Installation
=========
Change to directory DFTBopt and type 

./install.sh 

Currently, DFTBopt was tested for two compilers: GNU(g++) and INTEL(icpc).
You can choose the compiler by setting $CXX=“g++” or $CXX=“icpc” in the “install.sh” file.
To compile “erepopt”, a MPI library is also required. The install.sh script will try to get all the
required library and compile the code. You might have to adjust some flags in the makefile.
After the compilation, if you use bash, you can add following command to you “.bashrc” 
or run it before using DFTBopt,

source path-to-DFTBopt-directory/dftbopt_on.bash

Document
=============

There is a detail manual under doc/manual folder

Example
=============

For repopt, there are two examples rep1.in and rep2.in under the examples folder. 
It is straight forward to run these examples:
cd to the examples folder and type

repopt rep1.in >& rep1.out

repopt rep2.in >& rep2.out

For erepopt, the exmamples are under construction.
