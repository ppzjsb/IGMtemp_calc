
#######################################################################
#  Look at end of file for a brief guide to the compile-time options. #
#######################################################################

#--------------------------------------- Options
#OPTS +=-DPLAW_UVB
OPTS +=-DSECONDARY
#OPTS +=-DNO_HE2_HEAT

#--------------------------------------- Select system
#SYSTYPE="macbook13"
#SYSTYPE="linux"
SYSTYPE="local"

ifeq ($(SYSTYPE),"local")
CC = gcc		
OPTIMIZE  = -O0  
CVODEINCL = -I./sundials-2.7.0/include/
CVODELIB  = -L./sundials-2.7.0/lib/ -lsundials_cvode -lsundials_nvecserial
endif


ifeq ($(SYSTYPE),"macbook13")
CC = gcc		
OPTIMIZE  = -O0  
CVODEINCL = -I/Users/jamesbolton/FILES/Libraries/sundials/sundials-2.7.0_instdir/include/
CVODELIB  = -L/Users/jamesbolton/FILES/Libraries/sundials/sundials-2.7.0_instdir/lib/ -lsundials_cvode -lsundials_nvecserial
endif


ifeq ($(SYSTYPE),"linux")
CC = gcc		
OPTIMIZE  = -O3	
CVODEINCL = -I/home/ppzjsb/FILES/Libraries/sundials/sundials-2.7.0_instdir/include
CVODELIB  = -L/home/ppzjsb/FILES/Libraries/sundials/sundials-2.7.0_instdir/lib -lsundials_cvode -lsundials_nvecserial -Xlinker -R -Xlinker /home/ppzjsb/FILES/Libraries/sundials/sundials-2.7.0_instdir/lib
endif

CFLAGS = $(CVODEINCL) $(OPTIMIZE) -Wall $(OPTS)
LIBS   = $(CVODELIB) -lm 


EXEC = IGMtemp_calc
OBJS = main.o cooling.o utils.o uvb.o gausslegendre.o elec_interp.o
INCL = parameters.h proto.h global_vars.h elec_interp.h



$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS): $(INCL)

.PHONY: clean
clean:
	rm -f $(OBJS) $(EXEC)
	cd ./secondaries/; rm -f logxe_*.dat tablesize.dat

.PHONY: tidy
tidy:
	rm -f $(OBJS) $(EXEC) *~
	cd ./secondaries/; rm -f logxe_*.dat tablesize.dat *~


##############################################################################
#
# This code computes the ionisation state and temperature of a
# hydrogen and helium gas parcel of fixed comoving density.  Parts of
# the cooling.c routine are based on P-Gadget-3 (Springel et al. 2005,
# MNRAS, 364, 1105)
#	
# Options at compile-time: from the list below, activate/deactivate
# the options that apply to your run.  If you modify any of these
# options, make sure that you recompile the whole code by typing "make
# clean; make".
#	
# See also parameters.h for additional code settings
# 
# Last updated by James S. Bolton, 15/12/21
#
#
# Options:
#
#	- PLAW_UVB	Use a power-law UVB instead of a look up table based on a UV 
# 			background synthesis model.
#
#	- SECONDARY	Include secondary ionisations by fast photo-electrons. See
#			Furlanetto & Johnson-Stoever 2010, MNRAS, 404, 1869.
#
#	- NO_HE2_HEAT	Code testing option that turns off the He-II photo-heating 
#
##############################################################################
