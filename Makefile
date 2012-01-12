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
FFLAGS= -O -fopenmp -D$(ext) -DPPNAME
#-DPPNAME
#bugbuster -DPPNAME
#FFLAGS = -g -debug variable_locations -inline_debug_info -CB -check all 
LDADD = -lsrmodels
endif

INCLUDE = -I/usr/local/include/srmodels

INFOBJ = infprec.o cosmopar.o binfspline.o inftools.o infinout.o \
	 infmatter.o infdilaton.o infbgmodel.o infpotential.o infsigma.o \
	 infbounds.o infsric.o infbgfunc.o infbg.o infbgspline.o \
	 inftorad.o infpert.o infpowspline.o

FFLAGS +=  $(INCLUDE)

infbackmain.$(ext): $(INFOBJ) infbackmain.o
	$(FC) $(FFLAGS) $(INFOBJ) infbackmain.o -o $@ $(LDADD)

infpertmain.$(ext): $(INFOBJ) infpertmain.o
	$(FC) $(FFLAGS) $(INFOBJ) infpertmain.o -o $@ $(LDADD)

perttest.$(ext): $(INFOBJ) perttest.o
	$(FC) $(FFLAGS) $(INFOBJ) perttest.o -o $@ $(LDADD)

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90 

%.o: %.F90
	$(FC) $(FFLAGS) -c $*.F90 

clean:
	rm *.$(ext) *.o *.mod


