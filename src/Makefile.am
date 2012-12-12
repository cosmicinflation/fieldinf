SRC = infprec.f90 infunits.F90 binfspline.f90 infsolvers.f90 \
	infio.f90 infmatter.F90 infdilaton.f90 \
	infbgmodel.F90 infpotential.f90 infsigma.f90 \
	infbounds.f90 infsric.F90 infbgfunc.f90 infbg.f90 \
	infbgspline.f90 inftorad.f90 infpert.f90 infpowspline.f90

MOD = bspline.$(FC_MODEXT) infbg.$(FC_MODEXT) infbounds.$(FC_MODEXT) \
	infmatter.$(FC_MODEXT) infpowspline.$(FC_MODEXT) \
	infsric.$(FC_MODEXT) numeric.$(FC_MODEXT) infunits.$(FC_MODEXT) \
	infbgmodel.$(FC_MODEXT) infdilaton.$(FC_MODEXT) \
        infpert.$(FC_MODEXT) infprec.$(FC_MODEXT) infsolvers.$(FC_MODEXT) \
	infbgfunc.$(FC_MODEXT) infbgspline.$(FC_MODEXT) \
        infio.$(FC_MODEXT) infpotential.$(FC_MODEXT) \
	infsigma.$(FC_MODEXT) inftorad.$(FC_MODEXT)

check_PROGRAMS = fieldmain gurumain
fieldmain_SOURCES = $(SRC) fieldmain.f90
fieldmain_FCFLAGS = $(FCFLAGS)
fieldmain_LDFLAGS = $(LDFLAGS)
gurumain_SOURCES = $(SRC) gurumain.f90
gurumain_FCFLAGS = $(FCFLAGS)
gurumain_LDFLAGS = $(LDFLAGS)

lib_LTLIBRARIES = libfieldinf.la
libfieldinf_la_SOURCES = $(SRC)
libfieldinf_la_LDFLAGS = -version-info @version_info@
libfieldinf_la_includedir = $(includedir)/$(FIINCDIR)
libfieldinf_la_include_HEADERS = $(MOD)

man_MANS =

nodist_EXTRA_libfieldinf_la_SOURCES = donotexist.f90


clean-local: clean-modules clean-outfiles
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)
clean-outfiles:
	test -z "$(DATEXT)" || $(RM) *.$(DATEXT)