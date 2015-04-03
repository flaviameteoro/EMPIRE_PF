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















echo "Moving to empire base directory"
test cd ../


echo "Making empire codes"
test make > /dev/null

echo "Making minimal examples"
test make minimal > /dev/null

echo "Moving back to test directory"
test cd $dir

echo "Running comms initialise tests"
MPIRUN=/usr/bin/mpirun
MPIRUNOPTS="--output-filename out"


for mdls in {1..5}
do
    mdlprocs=$(($mdls * 5))
    for emps in $(eval echo {1..$mdls})
    do
	for i in $(eval echo {1..$mdls})
	do
	    echo "-np 5 ../bin/minimal_model_comms_v2 : " >> stuff
	done
	for j in $(eval echo {1..$emps})
	do
	    echo "-np 1 ../bin/minimal_empire_comms : " >> stuff
	done
	line=$(cat stuff | shuf | tr '\n' ' ' | rev | cut -c 4- | rev)
	test $MPIRUN $MPIRUNOPTS $line
	rm stuff
	outcheck
    done
done

echo "Running comms sending tests"
echo 2 > timesteps
for mdls in {1..5}
do
    mdlprocs=$(($mdls * 5))
    for emps in $(eval echo {1..$mdls})
    do
	for i in $(eval echo {1..$mdls})
	do
	    echo "-np 5 ../bin/minimal_model_v2 : " >> stuff
	done
	for j in $(eval echo {1..$emps})
	do
	    echo "-np 1 ../bin/minimal_empire : " >> stuff
	done
	line=$(cat stuff | shuf | tr '\n' ' ' | rev | cut -c 4- | rev)
	testp $MPIRUN $MPIRUNOPTS $line
	rm stuff
	outcheck
    done
done
rm timesteps


echo All tests of comms completed successfully.

