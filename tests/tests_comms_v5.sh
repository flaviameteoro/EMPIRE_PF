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










grep=/bin/grep
NC='\033[0m'
dir=`pwd`
MPIRUN=/usr/bin/mpirun
MPIRUNOPTS="--output-filename out"

commsfile=../comm_version.f90


version=$($grep comm_version= $commsfile | cut -f2 -d= | cut -f1 -d' ')
red EMPIRE VERSION detected as $version in $commsfile
if [[ "$version" = "5" ]]; then
    minmdl=../bin/minimal_model_v5
    minmdlcomm=../bin/minimal_model_comms_v5
    minmdlnum=5
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




testp $MPIRUN $MPIRUNOPTS -np 5 ../bin/minimal_model_comms_v5 : -np 1 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 10 ../bin/minimal_model_comms_v5 : -np 1 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 1 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 1 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 1 ../bin/minimal_empire_comms
rm out.*


testp $MPIRUN $MPIRUNOPTS -np 5 ../bin/minimal_model_comms_v5 : -np 4 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 10 ../bin/minimal_model_comms_v5 : -np 4 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 4 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 4 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 4 ../bin/minimal_empire_comms
rm out.*


testp $MPIRUN $MPIRUNOPTS -np 5 ../bin/minimal_model_comms_v5 : -np 2 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 5 ../bin/minimal_model_comms_v5 : -np 3 ../bin/minimal_empire_comms
rm out.*

testp $MPIRUN $MPIRUNOPTS -np 10 ../bin/minimal_model_comms_v5 : -np 2 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 10 ../bin/minimal_model_comms_v5 : -np 3 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 10 ../bin/minimal_model_comms_v5 : -np 5 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 10 ../bin/minimal_model_comms_v5 : -np 6 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 10 ../bin/minimal_model_comms_v5 : -np 7 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 10 ../bin/minimal_model_comms_v5 : -np 8 ../bin/minimal_empire_comms
rm out.*

testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 2 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 3 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 5 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 6 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 7 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 8 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 9 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 10 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 11 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 15 ../bin/minimal_model_comms_v5 : -np 12 ../bin/minimal_empire_comms
rm out.*

testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 2 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 3 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 5 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 6 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 7 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 8 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 9 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 10 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 11 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 12 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 13 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 14 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 15 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 20 ../bin/minimal_model_comms_v5 : -np 16 ../bin/minimal_empire_comms
rm out.*


testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 2 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 3 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 5 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 6 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 7 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 8 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 9 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 10 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 11 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 12 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 13 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 14 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 15 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 16 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 17 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 18 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 19 ../bin/minimal_empire_comms
rm out.*
testp $MPIRUN $MPIRUNOPTS -np 25 ../bin/minimal_model_comms_v5 : -np 20 ../bin/minimal_empire_comms
rm out.*






echo All tests of comms v5 completed successfully.

