#!/bin/bash
set -o nounset
NC='\033[0m'
test()
{
    blue $*
    $*
    if [[ $? != 0 ]]; then
	    echo "ERROR in that command. Return code $?. Stopping tests."
        if [[ -f model_specific.f90_backup ]]; then
            mv model_specific.f90_backup model_specific.f90
        fi
        if [[ -f ../model_specific.f90_backup ]]; then
            mv ../model_specific.f90_backup ../model_specific.f90
        fi
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
green()
{
    echo -e "\033[0;44m$* ${NC}"
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




testp cd ..
testp mkdir -p rundirectory
testp cp model_specific.f90 model_specific.f90_backup
testp cp examples/model_specific/identity/model_specific_identity.f90 model_specific.f90 
test make
test cd rundirectory


cat <<EOF > empire.nml
&pf_params
time_obs=4,
time_bwn_obs=50,
gen_data=.true.,
filter='DE',
init='N'
/

&comms_v4
nens=1
/

EOF

test $MPIRUN $MPIRUNOPTS -np 1 ../bin/empire
blue "Looking for observation files:"
ls -ltr obs*

if [[ ! -f obs_ts_0000050 ]]; then
    red There should be an observation file found. Stopping test.
    exit -1
fi

cat <<EOF > empire.nml
&pf_params
time_obs=4,
time_bwn_obs=50,
gen_data=.false.,
filter='DE',
init='N'
/

&comms_v4
nens=10
/

EOF
green Running 10 ensemble members on one process
time test $MPIRUN $MPIRUNOPTS -np 1 ../bin/empire

sed -i 's/nens=10/nens=100/g' empire.nml
green Running 100 ensemble members on one process
time test $MPIRUN $MPIRUNOPTS -np 1 ../bin/empire

sed -i 's/nens=100/nens=1000/g' empire.nml
green Running 1000 ensemble members on one process
time test $MPIRUN $MPIRUNOPTS -np 1 ../bin/empire

green Running 1000 ensemble members on 2 processes
time test $MPIRUN $MPIRUNOPTS -np 2 ../bin/empire

green Running 1000 ensemble members on 3 processes
time test $MPIRUN $MPIRUNOPTS -np 3 ../bin/empire

green Running 1000 ensemble members on 4 processes
time test $MPIRUN $MPIRUNOPTS -np 4 ../bin/empire

green Running 1000 ensemble members on 5 processes
time test $MPIRUN $MPIRUNOPTS -np 5 ../bin/empire

green Running 1000 ensemble members on 6 processes
time test $MPIRUN $MPIRUNOPTS -np 6 ../bin/empire

green Running 1000 ensemble members on 7 processes
time test $MPIRUN $MPIRUNOPTS -np 7 ../bin/empire

green Running 1000 ensemble members on 8 processes
time test $MPIRUN $MPIRUNOPTS -np 8 ../bin/empire

green Running 1000 ensemble members on 9 processes
time test $MPIRUN $MPIRUNOPTS -np 9 ../bin/empire

green Running 1000 ensemble members on 10 processes
time test $MPIRUN $MPIRUNOPTS -np 10 ../bin/empire

green Running 1000 ensemble members on 11 processes
time test $MPIRUN $MPIRUNOPTS -np 11 ../bin/empire

sed -i "s/filter='DE'/filter='SI'/g" empire.nml
purple Testing with the SIR filter

green Running 1000 ensemble members on one process
time test $MPIRUN $MPIRUNOPTS -np 1 ../bin/empire

green Running 1000 ensemble members on 2 processes
time test $MPIRUN $MPIRUNOPTS -np 2 ../bin/empire

green Running 1000 ensemble members on 3 processes
time test $MPIRUN $MPIRUNOPTS -np 3 ../bin/empire

green Running 1000 ensemble members on 4 processes
time test $MPIRUN $MPIRUNOPTS -np 4 ../bin/empire

green Running 1000 ensemble members on 5 processes
time test $MPIRUN $MPIRUNOPTS -np 5 ../bin/empire

green Running 1000 ensemble members on 6 processes
time test $MPIRUN $MPIRUNOPTS -np 6 ../bin/empire

green Running 1000 ensemble members on 7 processes
time test $MPIRUN $MPIRUNOPTS -np 7 ../bin/empire

green Running 1000 ensemble members on 8 processes
time test $MPIRUN $MPIRUNOPTS -np 8 ../bin/empire

green Running 1000 ensemble members on 9 processes
time test $MPIRUN $MPIRUNOPTS -np 9 ../bin/empire

green Running 1000 ensemble members on 10 processes
time test $MPIRUN $MPIRUNOPTS -np 10 ../bin/empire

green Running 1000 ensemble members on 11 processes
time test $MPIRUN $MPIRUNOPTS -np 11 ../bin/empire


if [[ -f ../model_specific.f90_backup ]]; then
    mv ../model_specific.f90_backup ../model_specific.f90
fi

purple All tests of comms_v4_full.sh completed successfully.
