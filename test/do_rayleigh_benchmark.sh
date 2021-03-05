#!/bin/bash
set -e
set -x

FLOTSAM=../bin/flotsam
OUTCODE="flotsam0522_rayleigh_od_"

WAVELENGTH=(543.8e-9 459.3e-9 368.3e-9 312.6e-9)
RAYLEIGH_OD=(0.1 0.2 0.5 1.0)

for INDEX in 0 1 2 3
do
    $FLOTSAM wavelength=${WAVELENGTH[$INDEX]} rayleigh_benchmark.cfg \
	> $OUTCODE${RAYLEIGH_OD[$INDEX]}_out.dat
done

echo "TOA FLUXES:"
for INDEX in 0 1 2 3
do
    echo ${RAYLEIGH_OD[$INDEX]} $(sed -n 2p $OUTCODE${RAYLEIGH_OD[$INDEX]}_out.dat | awk '{print $10}')
done
