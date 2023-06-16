
---------------------------------------------------------------

Brief notes on use of the non-equilibrium photoionisation code used in
Bolton et al. (2022), MNRAS, 513, 864
Bolton, Caputo, Liu & Viel (2022), PRL, 129, 211102

---------------------------------------------------------------

RUNNING THE CODE:

Requirements are a C compiler and the numerical library (CVODE) used
for solving the stiff differential equations for the thermochemistry.
For this, the SUNDIALS v2.7.0 package must be installed.  Follow the
instructions at this link to install the library.

https://computing.llnl.gov/projects/sundials/sundials-software

To compile the code, at the command line type:
make

To run the code, at the command line type:
./IGMtemp_calc

A binary output file will then appear in the output/ directory.

Parameters that can be varied in the code are contained in the file
parameters.h.

Compile time options are included in the Makefile, along with descriptions.

If changing the compile time options one should generate a new binary.
To do this, at the command line type:
make clean
make


ADDITIONAL FILES:

uvb_models/ contains a selection of ASCII tables with the
photoionisation and photoheating rates predicted by a variety of
published UVB synthesis models.  These can be selected within
parameters.h.  The reference papers are:

Haardt & Madau, 2012, ApJ, 746, 125
Khaire & Srianand, 2019, MNRAS, 484, 4174
Puchwein et al. 2019, MNRAS, 485, 47
Faucher-Giguere 2020, MNRAS, 493, 1614

python_routines/ contains a simple python code that reads in the
output from the C code and plots a selection of data.

idl_routines/ contains a simple IDL code that does the same thing as
the python routine. 

secondaries/ contains binary look-up tables for the implementation of
secondary ionisations by fast electrons.  The reference paper is:

Furlanetto & Johnson-Stoever, 2010, MNRAS, 404, 1869

---------------------------------------------------------------

