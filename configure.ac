AC_INIT([fieldinf], [1.0.12], [chris@ringeval.com])

AC_PROG_FC
AC_PROG_CPP
AC_PROG_LIBTOOL

version_info="0:1:0"
AC_SUBST(version_info)


FCFLAGS+=" -fopenmp"

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

AX_F90_LIBRARY_SETUP([aspic],[],[cosmopar.$FC_MODEXT],[libaspic.so],
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

AM_INIT_AUTOMAKE

AC_CONFIG_FILES([Makefile
                 src/Makefile])                 

AC_OUTPUT