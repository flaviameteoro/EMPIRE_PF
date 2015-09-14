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






NC='\033[0m'
dir=`pwd`
grep=/bin/grep
commsfile=../src/utils/comms.f90


for version in {1..3}; do

    echo "Moving to empire base directory"
    test cd ../

    echo "Changing to empire version $version"
    test make v$version

    echo "Cleaning empire build"
    test make clean

    echo "Compiling all codes"
    test make all

    red Compilation of empire codes for version $version successful

    testp cd tests

    echo "running test of the comms then"
    test test_comms.sh
done

purple All tests of comms completed successfully for all versions

