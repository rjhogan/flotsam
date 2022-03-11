#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])

# Remember to update the version number for a new release
AC_INIT([FLOTSAM], [0.5.23], [r.j.hogan@ecmwf.int],[flotsam])

# Most code is C++ but there is some C and some Fortran
AC_LANG([C++])

AC_CONFIG_SRCDIR([include/flotsam.h])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Wno-error foreign])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CXX
#AC_PROG_F77
#AC_PROG_FC
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
#AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Check for OpenMP
AC_OPENMP
AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")

# Header directories needed during the compilation
AC_SUBST(AM_CPPFLAGS,"-I../include")

# Check for NetCDF library: REQUIRED
#AC_SEARCH_LIBS([nc_create], [netcdf], , AC_MSG_ERROR([Unable to find NetCDF library]))

# Check for Adept library: REQUIRED
LIBS="$LIBS -ladept"
AC_MSG_CHECKING([for Adept >= 2.0.4: including adept_arrays.h and linking via -ladept])
AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <adept_arrays.h>
#include <string>
#if ADEPT_VERSION < 20004
#error "Adept version >= 2.0.4 required"
#endif],[std::string test = adept::compiler_version()])],AC_MSG_RESULT([yes]),AC_MSG_RESULT([no])
AC_MSG_ERROR([Unable to find Adept library version >= 2.0.4]))

# Checks for standard header files
#AC_CHECK_HEADERS([fenv.h malloc.h memory.h stdint.h stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics
AC_TYPE_SIZE_T

# Due to autotools stupidity, this breaks the build:
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([floor gethostname gettimeofday memset pow sqrt strcasecmp strchr strdup strncasecmp strstr strtol])

# Makefiles to create
AC_CONFIG_FILES([Makefile
		 include/Makefile
                 src/Makefile
		 bin/Makefile
		 test/Makefile])

# Provide a summary to the user
AC_MSG_NOTICE([********************* Summary **************************************])
AC_MSG_NOTICE([  CXX = $CXX ])
AC_MSG_NOTICE([  CPPFLAGS = $CPPFLAGS])
AC_MSG_NOTICE([  CXXFLAGS = $CXXFLAGS])
AC_MSG_NOTICE([  LDFLAGS =  $LDFLAGS])
AC_MSG_NOTICE([  LIBS = $LIBS])
AC_MSG_NOTICE([Typing "make; make install" will install the code and auxiliary files])
AC_MSG_NOTICE([********************************************************************])

AC_OUTPUT
