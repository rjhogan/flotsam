FLOTSAM=../bin/flotsam
VERS=0521a_fit
CFGFILE=one_profile.cfg
DIRECTORY=$HOME/work/monte_carlo/data
PFCODE=du_1024nm
#PFCODE=is_0470nm
PROPCFGFILE=CAMS_benchmark/properties_${PFCODE}.cfg
PFFILE=CAMS_benchmark/phase_function_${PFCODE}.txt

PFCODE=rayleigh
PFCODE=isotropic
#PFCODE=droplet_re10um_550nm
PROPCFGFILE=
PFFILE=phase_function_${PFCODE}.txt

set -x

for SZA in 0 20 40 60 80
do
    for OD in 0.5 2 10
    do
	FILEBASE=$(echo internals_${VERS}_${PFCODE}_sza${SZA}_od${OD} | tr '.' '_')
	$FLOTSAM optical_depth=$OD sza=$SZA mode=internals \
	    phase_function_file=${PFFILE} \
	    $PROPCFGFILE $CFGFILE > $DIRECTORY/$FILEBASE.m
    done
done

$FLOTSAM phase_function_file=${PFFILE} droplet10micron.cfg > $SCRATCH/flotsam/flotsam${VERS}_nogas_0.00_${PFCODE}_out.dat
