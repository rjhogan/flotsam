#!/bin/bash

# This script calls the ./configure script with settings that should
# provide the best performance with the GNU compiler. Any additional
# arguments are passed to the configure script.

CXX="g++"
CXXFLAGS="-Wall -g1 -O3 -march=native"
CPPFLAGS="-DADEPT_NO_ALIAS_CHECKING -DADEPT_NO_DIMENSION_CHECKING -DADEPT_STACK_THREAD_UNSAFE" 

if [ ! -z "$ADEPT_PREFIX" ]
then
    # If the ADEPT_PREFIX environment variable is set then ensure
    # Adept include files and shared library can be located
    if [ -d "$ADEPT_PREFIX/lib64" ]
    then
	LDFLAGS="-L$ADEPT_PREFIX/lib64 -Wl,-rpath,$ADEPT_PREFIX/lib64"
    else
	LDFLAGS="-L$ADEPT_PREFIX/lib -Wl,-rpath,$ADEPT_PREFIX/lib"
    fi
    CPPFLAGS+=" -I$ADEPT_PREFIX/include"
    set -x
    ./configure "CXX=$CXX" "CXXFLAGS=$CXXFLAGS" "CPPFLAGS=$CPPFLAGS" "LDFLAGS=$LDFLAGS" $@

else
    # Adept is assumed to be in the system search path
    set -x
    ./configure "CXX=$CXX" "CPPFLAGS=$CPPFLAGS" "CXXFLAGS=$CXXFLAGS" $@

fi



