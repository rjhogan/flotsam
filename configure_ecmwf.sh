#!/bin/bash

# This script calls the ./configure script with settings appropriate for ECMWF

set -x

# Set install location
#INSTALLDIR=/var/tmp/$USER
INSTALLDIR=$HOME/apps/flotsam-0.5.22

# Select version and location of Adept
ADEPT_VER=adept-2.1
ADEPTHOME=/home/rd/parr
ADEPT_PREFIX=$ADEPTHOME/apps/$ADEPT_VER
ADEPT_LDFLAGS="-L$ADEPT_PREFIX/lib64 -Wl,-rpath,$ADEPT_PREFIX/lib64"

# Prepare variables to pass to configure script
CPPFLAGS="-I$ADEPT_PREFIX/include"
#CPPFLAGS="-I/home/rd/parr/git/adept2/include"
#CPPFLAGS="-I$ADEPT_PREFIX/include -DADEPT_BOUNDS_CHECKING"

LDFLAGS="$ADEPT_LDFLAGS" 

CXXFLAGS="-Wall -g1 -O3 -march=native -std=c++11"
#CXXFLAGS="-Wall -g1 -O3 -march=native -std=c++11 -DADEPT_FAST" 
#CXXFLAGS="-Wall -g -O3 -std=c++11 -DADEPT_BOUNDS_CHECKING"

# Call configure script
./configure "CXXFLAGS=$CXXFLAGS" "CPPFLAGS=$CPPFLAGS" "LDFLAGS=$LDFLAGS" \
    --prefix "$INSTALLDIR" $@
