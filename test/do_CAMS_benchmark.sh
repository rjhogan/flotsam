#!/bin/bash

# This script runs the CAMS benchmark cases for various wavelengths,
# surfaces, aerosol types and viewing geometries

# Exit on failure
set -e

# All 49 layers
MODE=full

# Coarsen the profile to 12 layers
#MODE=coarse

# Simplify the calculation: (1) no molecular absorption, (2) Rayleigh
# and aerosol scattering properties are assumed to be perfectly mixed
# throughout the profile, (3) no aerosol optical depth scaling is
# performed to account for the requested optical depth being at 550nm
# which may be different from the actual wavelength.
#MODE=simple

# Use the actual pressure profile but turn off gas absorption or gas
# absorption and scattering (note that FLOTSAM will fail in the latter
# case if the AOD=0 because then the atmosphere is a complete
# vacuum). In both these cases the optical_depth_scaling=1 for all
# wavelengths.
#MODE=nogasabs
#MODE=nogas

EXTRAS="mode=direct"
#EXTRAS="mode=jacobian"

# Input data to be found here
CAMSDIR="CAMS_benchmark"

# Output directory
PREFIX="."
#PREFIX="/scratch/rd/parr/flotsam"

# Output files prefixed by this
OUTCODE="flotsam0522_$MODE"

# Executable
FLOTSAM=../bin/flotsam

# Configuration file matching the viewing geometry required by the
# CAMS benchmark intercomparison
MASTER_FILE=${CAMSDIR}/CAMS_benchmark.cfg

# Alternative configuration file with more instrument elevations and
# azimuths, for plotting the results
#MASTER_FILE=${CAMSDIR}/CAMS_benchmark_plot.cfg

# The profile can be coarsened by specifying a list of integers
# representing the numbers of layers to group, starting from the top,
# with any remaining layers left ungrouped
if [ "$MODE" = "coarse" ]
then
    COARSENING='coarsening=38'
else
    COARSENING=
fi

if [ "$MODE" = nogas ]
then
    EXTRAS="$EXTRAS no_rayleigh=1"
fi

if [ "$MODE" = nogasabs -o "$MODE" = nogas ]
then
    EXTRAS="$EXTRAS optical_depth_scaling=1"
fi

# In the simple configuration we override the number of layers, and
# the optical depth scaling of 1 ensures (3) is satisfied above
SIMPLE_CONFIG="layer_count=49 optical_depth_scaling=1 surface_pressure=101300.0"

# Loop over wavelengths in nm
for WAVELENGTH in 0470 0550 0660 0865 1024
do
    if [ "$MODE" = nogasabs -o "$MODE" = nogas ]
    then
	MOLECULAR_FILE=
    else
	MOLECULAR_FILE="${CAMSDIR}/CAMS_benchmark_molecular_optical_depth_${WAVELENGTH}.cfg"
    fi

    # Loop over aerosol types dust, sea-salt, industrial absorbing and
    # industrial scattering
    for TYPE in du is ss ia
    do
	# Loop over surface types: Lambertian with three albedo
	# values, and Cox-Munk
	for ALBEDO in 0.00 0.05 0.10 coxm
	do
	    if [ "$ALBEDO" = "coxm" ]
	    then
		SURFACE="surface=cox_munk"
	    else
		SURFACE="albedo=$ALBEDO"
	    fi
	    OUTFILE="${PREFIX}/${OUTCODE}_${ALBEDO}_${TYPE}_${WAVELENGTH}nm_out.dat"

	    if [ "$MODE" = simple ]
	    then
		# In the simple case the
		# CAMS_benchmark_molecular_optical_depth_*.cfg file is
		# not used
		set -x
		$FLOTSAM $EXTRAS $SURFACE \
		    phase_function_file=${CAMSDIR}/phase_function_${TYPE}_${WAVELENGTH}nm.txt \
		    wavelength=${WAVELENGTH}e-9 ${SIMPLE_CONFIG}\
		    ${MASTER_FILE} ${CAMSDIR}/properties_${TYPE}_${WAVELENGTH}nm.cfg \
		    > ${OUTFILE}
		set +x
		echo "WRITTEN ${OUTFILE}"
	    else
		set -x
		$FLOTSAM $EXTRAS $SURFACE $COARSENING wavelength=${WAVELENGTH}e-9 \
		    phase_function_file=${CAMSDIR}/phase_function_${TYPE}_${WAVELENGTH}nm.txt \
		    ${MASTER_FILE} \
		    ${CAMSDIR}/CAMS_benchmark_aerosol_optical_depth.cfg \
		    ${MOLECULAR_FILE} \
		    ${CAMSDIR}/properties_${TYPE}_${WAVELENGTH}nm.cfg \
		    > ${OUTFILE}
		set +x
		echo "WRITTEN ${OUTFILE}"
	    fi
	done
    done
done
