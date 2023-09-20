#!/bin/bash

# This script calls the ./configure script with settings appropriate for ECMWF

set -x

# Set install location
#INSTALLDIR=/var/tmp/$USER
INSTALLDIR=$HOME/apps/flotsam-0.5.23

# Select version and location of Adept
#ADEPT_VER=adept-2.1-sparse-nothreads
ADEPT_VER=adept-2.1.2
ADEPTHOME=/home/parr
ADEPT_PREFIX=$ADEPTHOME/apps/$ADEPT_VER
ADEPT_LDFLAGS="-L$ADEPT_PREFIX/lib -Wl,-rpath,$ADEPT_PREFIX/lib"

# Prepare variables to pass to configure script
CPPFLAGS="-I$ADEPT_PREFIX/include -DADEPT_STORAGE_THREAD_SAFE"
#CPPFLAGS="-I/home/rd/parr/git/adept2/include"
#CPPFLAGS="-I$ADEPT_PREFIX/include -DADEPT_BOUNDS_CHECKING"

LDFLAGS="$ADEPT_LDFLAGS -L$HOME/apps/openblas-0.3.17-zen/lib" 

# Optimized
CXXFLAGS="-Wall -g1 -O2 -march=native -std=c++11"
CXXFLAGS="-Wall -g1 -O3 -march=native -std=c++11"
#CXXFLAGS="-Wall -g1 -O2 -std=c++11"
#CXXFLAGS="-Wall -g1 -O2 -std=c++11 -pg"
# Fast but possibly unsafe
#CXXFLAGS="-Wall -g1 -O3 -march=native -std=c++11 -DADEPT_FAST"
# Debug version
#CXXFLAGS="-Wall -g -Og -std=c++11 -DADEPT_BOUNDS_CHECKING"
# Debug unoptimized version
#CXXFLAGS="-Wall -g -O0 -std=c++11 -DADEPT_BOUNDS_CHECKING"

# Call configure script
./configure "CXXFLAGS=$CXXFLAGS" "CPPFLAGS=$CPPFLAGS" "LDFLAGS=$LDFLAGS" \
	    --prefix "$INSTALLDIR"
#CXX=g++-12 CC=gcc-12$@
