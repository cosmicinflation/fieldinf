AC_INIT([fieldinf], [1.0.8], [chris@ringeval.com])

AC_PROG_FC
AC_PROG_CPP
AC_PROG_LIBTOOL

version_info="0:0:0"
AC_SUBST(version_info)


FCFLAGS=" -O -fopenmp -DPPNAME"

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

AX_F90_LIBRARY_SETUP([aspic],[],[cosmopar.$FC_MODEXT],[libaspic*],
[
	use cosmopar
	write(*,*)'kstar= ',kstar
])

AM_CONDITIONAL([USE_APIC], [test x$F90_LIBS_aspic != x])
AM_COND_IF (USE_ASPIC,
[
	LDFLAGS="$(echo "$LDFLAGS" | sed 's/-L /-L/g')"

],
[
	FCFLAGS="$(echo "FCFLAGS" | sed 's/-DPPNAME//g')"
	FCFLAGS+=" -DNOASPIC"
	AM_MESG_WARN([libaspic not found!])
	AM_MESG_WARN([you can specify a path using: --with-aspic=])
	AM_MESG_WARN([libaspic support DISABLED])
])	

AM_INIT_AUTOMAKE

AC_CONFIG_FILES([Makefile
                 src/Makefile])                 

AC_OUTPUT