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







here=$(pwd)
top=$(echo $here | rev | cut -f1 -d/ | rev)
if [[ "$top" != "tests" ]]; then
    test cd tests
fi





grep=/bin/grep

dir=`pwd`
MPIRUN=/usr/bin/mpirun
MPIRUNOPTS="--output-filename out"
commsfile=../src/utils/comms.f90


version=$($grep empire_version= $commsfile | cut -f2 -d=)
red EMPIRE VERSION detected as $version in $commsfile
if [[ "$version" = "1" ]]; then
    minmdl=$(readlink -e ../bin/lorenz96)
elif [[ "$version" = "2" ]]; then
    minmdl=$(readlink -e ../bin/lorenz96_v2)
else
    red ERROR: version not supported by these tests.
    exit -1
fi
empire=$(readlink -e ../bin/empire)
echo $empire
echo $minmdl





echo "Moving to empire base directory"
test cd ../


echo "Making empire codes"
test make > /dev/null

echo "Making lorenz96 model"
test make lorenz96 > /dev/null

echo "Moving back to temporary test directory"
tempdir=$(mktemp -d)
test cd $tempdir



echo "running gen obs"
line="-np 1 $minmdl : -np 1 $empire"
cat <<EOF > empire.nml
&pf_params
time_obs=1,
time_bwn_obs=4,
gen_data=.true.,
init='N'
/
EOF
test $MPIRUN $MPIRUNOPTS $line


cat <<EOF > empire.nml
&pf_params
time_obs=1,
time_bwn_obs=4,
filter='LD',
init='N'
/
EOF
echo "running LETKF"
line="-np 8 $minmdl : -np 2 $empire"
test $MPIRUN $MPIRUNOPTS $line

echo All tests of lorenz96 completed successfully.

