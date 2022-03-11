
---------------------------------------------------------------

Brief notes on use of the non-equilibrium photoionisation code used in
Bolton et al. (2021), MNRAS submitted, arXiv:2111.09600

James S. Bolton,
Nottingham, 11/03/22

---------------------------------------------------------------

Compile time options are included in the Makefile.  Requirements are a
C compiler and link to the numerical library used for solving the
stiff differential equations for the ionisation balance (SUNDIALS
v2.7.0, my local build has been provided in this bundle).

https://computing.llnl.gov/projects/sundials/sundials-software

Warning: it may be the case that you need to rebuild your own version
of the SUNDIALS/CVODE library to run the code.  Follow the instructions at
the above website.

To compile the code, at the command line type:
make

To run the code, at the command line type:
./IGMtemp_calc

A binary output file will then appear in the output/ directory.

Parameters that can be varied in the code are contained in the file
parameters.h.  

uvb_models/ contains a selection of ASCII tables with the
photoionisation and photoheating rates predicted by a variety of
published UVB synthesis models.  These can be selected within
parameters.h

python_routines/ contains a simple python code that reads in the
output from the C code and plots the data within.

secondaries/ contains binary look-up tables for the implementation of
secondary ionisations by fast electrons.  Note that secondary
ionisations have a very small effect on the gas temperature at z~0.

---------------------------------------------------------------

