# >>> DESIGNED FOR GMAKE <<<

FC=f95

ext=$(shell uname | cut -c1-3)


ifeq ($(ext),IRI)
FFLAGS= -O -64 -OPT:Olimit=0 -D$(ext)
endif

ifeq ($(ext),OSF)
FFLAGS= -O  -D$(ext)
endif

ifeq ($(ext),Lin)
FC=gfortran
FFLAGS= -O -fopenmp -D$(ext) -DPP12
#bugbuster
#FFLAGS = -g -debug variable_locations -inline_debug_info -CB -check all 
endif

INFOBJ = infprec.o cosmopar.o binfspline.o inftools.o hyper2F1.o specialinf.o infinout.o infbgmodel.o \
	infsrmodel.o infbg.o infbgspline.o inftorad.o infpert.o infpowspline.o


infbackmain.$(ext): $(INFOBJ) infbackmain.o
	$(FC) $(FFLAGS) $(INFOBJ) infbackmain.o -o $@

infpertmain.$(ext): $(INFOBJ) infpertmain.o
	$(FC) $(FFLAGS) $(INFOBJ) infpertmain.o -o $@

perttest.$(ext): $(INFOBJ) perttest.o
	$(FC) $(FFLAGS) $(INFOBJ) perttest.o -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90

%.o: %.F90
	$(FC) $(FFLAGS) -c $*.F90

clean:
	rm *.$(ext) *.o *.mod


