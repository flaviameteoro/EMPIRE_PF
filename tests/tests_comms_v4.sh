#!/bin/bash
set -o nounset

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
        if [[ -f empire.nml ]]; then
            rm empire.nml
        fi
	exit -1
    fi
}
testr()
{
    red The following command should break and return an error:
    red $*
    $*
    if [[ $? == 0 ]]; then
	    echo "ERROR in that command. Return code $?. Stopping tests."
        if [[ -f empire.nml ]]; then
            rm empire.nml
        fi
	    exit -1
    else
        red "and it did"
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

setens()
{
cat <<EOF>empire.nml
&comms_v4
nens=$1
/

EOF
cat empire.nml
}








grep=/bin/grep
NC='\033[0m'
dir=`pwd`
MPIRUN=/usr/bin/mpirun
MPIRUNOPTS="--output-filename out"

commsfile=../src/comm_version.f90


version=$($grep comm_version= $commsfile | cut -f2 -d= | cut -f1 -d' ')
red EMPIRE VERSION detected as $version in $commsfile
if [[ "$version" = "4" ]]; then
    blue Doing the tests for EMPIRE comms version 4
else
    red ERROR: version not supported by these tests.
    exit -1
fi


echo "Moving to empire base directory"
test cd ../


echo "Making empire codes"
test make > /dev/null

echo "Making minimal examples"
test make minimal > /dev/null

echo "Moving back to test directory"
test cd tests

echo "Running comms initialise tests"


setens 1
testp $MPIRUN $MPIRUNOPTS -np 1 ../bin/minimal_empire_comms
setens 2
testp $MPIRUN $MPIRUNOPTS -np 1 ../bin/minimal_empire_comms
setens 3
testp $MPIRUN $MPIRUNOPTS -np 1 ../bin/minimal_empire_comms
setens 4
testp $MPIRUN $MPIRUNOPTS -np 1 ../bin/minimal_empire_comms
setens 5
testp $MPIRUN $MPIRUNOPTS -np 1 ../bin/minimal_empire_comms

setens 2
testp $MPIRUN $MPIRUNOPTS -np 1 ../bin/minimal_empire_comms
setens 2
testp $MPIRUN $MPIRUNOPTS -np 2 ../bin/minimal_empire_comms
setens 3
testp $MPIRUN $MPIRUNOPTS -np 1 ../bin/minimal_empire_comms
setens 3
testp $MPIRUN $MPIRUNOPTS -np 2 ../bin/minimal_empire_comms
setens 3
testp $MPIRUN $MPIRUNOPTS -np 3 ../bin/minimal_empire_comms


setens 1
testr $MPIRUN $MPIRUNOPTS -np 2 ../bin/minimal_empire_comms


if [[ -f empire.nml ]]; then
    rm empire.nml
    rm out.1.*
fi
echo All tests of comms v4 completed successfully.

