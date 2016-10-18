#!/bin/bash
set -o nounset

NC='\033[0m'
test()
{
    blue $*
    $*
    if [[ $? != 0 ]]; then
	echo "ERROR in that command. Return code $?. Stopping tests."
	exit -1
    fi
}
testp()
{
    purple $*
    $*
    if [[ $? != 0 ]]; then
	echo "ERROR in that command. Return code $?. Stopping tests."
	mv $empdir/$ms_backup $empdir/model_specific.f90
	exit -1
    fi
}
blue()
{
    echo -e "\033[0;34m$* ${NC}"
}
white()
{
    echo -e "\033[1;37m$* ${NC}"
}
red()
{
    echo -e "\033[0;31m$* ${NC}"
}
purple()
{
    echo -e "\033[1;35m$* ${NC}"
}

outcheck()
{
    $grep -rni warning out*
    if [[ $? = 0 ]]; then
	echo 'WARNINGS FOUND IN OUTPUT FILES. STOPPING'
	exit -1
    fi
    $grep -rni error out*
    if [[ $? = 0 ]]; then
	echo 'ERRORS FOUND IN OUTPUT FILES. STOPPING'
	exit -1
    fi
    $grep -rni success > /dev/null
    if [[ $? != 0 ]]; then
	echo 'No successes found in output files'
    fi
    red Success
    rm out*
}




red TESTING 3DVar...

purple 'B = diag(2,3,4), R = I, H = I'

here=$(pwd)
top=$(echo $here | rev | cut -f1 -d/ | rev)
if [[ "$top" != "tests" ]]; then
    test cd tests
fi

blue Compiling prerequisites
cd ..
make EMPIRE models > /dev/null 2> /dev/null
cd -


grep=/bin/grep

dir=`pwd`
MPIRUN=/usr/bin/mpirun
MPIRUNOPTS="--output-filename out"
commsfile=../comm_version.f90
ms_backup=model_specific.f90_real




version=$($grep comm_version= $commsfile | cut -f2 -d=)
red EMPIRE VERSION detected as $version in $commsfile
if [[ "$version" = "1" ]]; then
    minmdl=$(readlink -e ../bin/linear_empire_vader)
elif [[ "$version" = "2" ]]; then
    minmdl=$(readlink -e ../bin/linear_empire_vader_v2)
else
    red ERROR: version not supported by these tests.
    exit -1
fi
empire=$(readlink -e ../bin/empire)

if [[ "$minmdl" = "" ]]; then
    echo "for some reason minmdl was not found. Stopping tests early"
    exit -2
fi

echo $empire
echo $minmdl





echo "Moving to empire base directory"
test cd ../



empdir=$(pwd)
if [ -f $ms_backup ]; then
    echo "file $ms_backup exists"
    echo "will not continue unless this is removed"
    exit -1
fi

echo "backing up model_specific.f90"
cp $empdir/model_specific.f90 $empdir/$ms_backup


white "copying the lorenz96 model_specfic file into place"
testp cp examples/linear/model_specific_linear_3dvar.f90 model_specific.f90

echo "Making empire codes"
testp make > /dev/null

echo "Making linear model"
testp make linear > /dev/null

echo "Moving back to temporary test directory"
tempdir=$(mktemp -d)
testp cd $tempdir



echo "running gen obs"
line="-np 1 $minmdl : -np 1 $empire"
cat <<EOF > empire.nml
&pf_params
time_obs=1,
time_bwn_obs=1,
gen_data=.true.,
init='N'
/

&empire_output
basename='truth'
/

EOF
MPIRUNOPTS="--output-filename truth"
testp $MPIRUN $MPIRUNOPTS $line


cat <<EOF > empire.nml
&pf_params
time_obs=1,
time_bwn_obs=1,
filter='3D',
init='N'
/

&var_params
opt_method='cg',
cg_method=2,
cg_eps=1.0d-5,
total_timesteps=1,
/

&empire_output
basename='cg'
/

EOF
echo "running 3DVar using nonlinear CG"
MPIRUNOPTS="--output-filename cg"
line="-np 4 $minmdl : -np 1 $empire"
testp $MPIRUN $MPIRUNOPTS $line




cat <<EOF > empire.nml
&pf_params
time_obs=1,
time_bwn_obs=1,
filter='3D',
init='N'
/

&var_params
opt_method='lbfgs',
total_timesteps=1,
/

&empire_output
basename='bfgs'
/

EOF
echo "running 3DVar using L-BFGS"
MPIRUNOPTS="--output-filename lbfgs"
line="-np 4 $minmdl : -np 1 $empire"
testp $MPIRUN $MPIRUNOPTS $line



echo All tests of 3DVar completed successfully.
mv $empdir/$ms_backup $empdir/model_specific.f90
