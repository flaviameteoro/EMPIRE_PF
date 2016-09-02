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
commsfile=../src/comm_version.f90




echo "Moving to empire base directory"
test cd $(git rev-parse --show-toplevel)

for version in {5..1..1}; do
    echo "Changing to empire version $version"
    test make v$version

    echo "Cleaning empire build"
    test make clean

    echo "Compiling all codes"
    test make all

    red Compilation of empire codes for version $version successful
done

purple All compilation of codes completed successfully.

