#!/bin/bash
set -e
set -x

FLOTSAM=../bin/flotsam
OUTCODE="flotsam0522_iso"

ISO_OD="0.001 0.01 0.1 0.2 0.5 1.0"

$FLOTSAM wavelength=500e-9 "optical_depth=$ISO_OD" single_scattering_albedo=1 no_rayleigh=1 phase_function_file=phase_function_isotropic.txt rayleigh_benchmark.cfg \
	> ${OUTCODE}_out.dat
