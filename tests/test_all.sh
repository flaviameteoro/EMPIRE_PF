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







here=$(pwd)
top=$(echo $here | rev | cut -f1 -d/ | rev)
if [[ "$top" != "tests" ]]; then
    test cd tests
fi


#compilation tests:
test test_compile.sh
test test_examples_ms.sh
test test_comms_all_versions.sh


cd ..;make v1 clean;cd tests

test test_l96.sh
test test_instructions_l63.sh
test test_instructions_l96.sh

red All tests done. Phew.
