AC_INIT([fieldinf],[1.2.0],[fieldinf@ringeval.info],[],[http://cp3.irmp.ucl.ac.be/~ringeval/fieldinf.html])

AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_FC
AC_PROG_CPP
AC_PROG_LIBTOOL

version_info="3:0:0"
AC_SUBST(version_info)


FCFLAGS+=" -fopenmp -DPPNAME"

FIINCDIR="fieldinf"
AC_SUBST([FIINCDIR])

AX_F90_MODULE_EXTENSION
if test x$ax_cv_f90_modext = xunknown ; then
   AC_MSG_ERROR([unable to find f90 modules extension])
fi

FC_MODEXT=$ax_cv_f90_modext
AC_SUBST([FC_MODEXT]) 

DATEXT="dat"
AC_SUBST([DATEXT]) 

AC_FC_SRCEXT(F90)

AX_F90_LIBRARY_SETUP([aspic],[aspic.h],[],[libaspic.so],
[
	use cosmopar
	write(*,*)'kstar= ',kstar
])

AM_CONDITIONAL([USE_ASPIC],[test x$F90_LIBS_aspic != x])
AM_COND_IF([USE_ASPIC],
[
	LDFLAGS="$(echo "$LDFLAGS" | sed 's/-L /-L/g')"
],
[
	FCFLAGS="$(echo "$FCFLAGS" | sed 's/-DPPNAME//g')"
	FCFLAGS+=" -DNOASPIC"
	AC_MSG_WARN([libaspic not found!])
	AC_MSG_WARN([you can specify a path using: --with-aspic=])
	AC_MSG_WARN([libaspic support DISABLED])
])	

AC_CONFIG_TESTDIR([tests])

AC_CONFIG_FILES([Makefile
                 tests/Makefile
                 src/Makefile])

AC_OUTPUT