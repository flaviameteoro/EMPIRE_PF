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
	mv $ms_backup model_specific.f90
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

dir=`pwd`
MPIRUN=/usr/bin/mpirun
MPIRUNOPTS="--output-filename out"
commsfile=../src/utils/comms.f90
ms_backup=model_specific.f90_real


echo Test of build with model specific files................



#directory check:
if [[ ! -f model_specific.f90 ]] ; then
    #echo your mum
#else
   if [[ -f ../model_specific.f90 ]] ; then
       test cd ../
   else
       red "ERROR: cannot find appropriate test directory"
   fi
fi









if [[ -f $ms_backup ]];then
    red $ms_backup exists
    red Will not continue unless this file is moved
    red ERROR
    exit -1
else
    blue "Copying model_specific.f90 to $ms_backup"
    test cp model_specific.f90 $ms_backup
fi




for ms in examples/*/model_specific*.f90
do
    echo $ms
    testp cp $ms model_specific.f90
    testp make clean
    testp make > /dev/null
    echo EMPIRE build successful for $ms
done
purple mv $ms_backup model_specific.f90
mv $ms_backup model_specific.f90

echo All tests of examples completed successfully.
exit 0
