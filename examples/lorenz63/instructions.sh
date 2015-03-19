#!/bin/bash
#instructions for testing/an example for 4dEnVar

######################bash preliminaries=============
set -o nounset
set -o errexit
comment()
{
echo -e "${red}$* ${NC}"
}

command()
{
echo -e "${blue}$* ${NC}"
}

red='\033[0;31m'
blue='\033[0;34m'
NC='\033[0m'
dir=`pwd`
####################################################

F90=mpif90
mdlex=lorenz63_empire_vader
empex=EMPIRE_4DENVAR

comment "#build the model:"
rm -rf $mdlex
command $F90 -o $mdlex Lorenz63_empire_vader.f90
$F90 -o $mdlex Lorenz63_empire_vader.f90
echo '#Model built successfully'

#move the model to a "run folder"
mkdir -p run/
mv $mdlex run/


comment '#put an appropriate model specific file in the right place:'
command cp ../../model_specific.f90 model_specific_backup.f90
cp ../../model_specific.f90 model_specific_backup.f90
command cp model_specific_l63.f90 ../../model_specific.f90
cp model_specific_l63.f90 ../../model_specific.f90


comment '#move to the root folder:'
command cd ../../
cd ../../

comment '#Build the full empire code:'
make clean
make > /dev/null 2>&1
command make
echo 'Built empire filter codes successfully'

comment '#copy the empire code to the run directory:'
command cp bin/empire $dir/run
cp bin/empire $dir/run/

comment '#move back to this example directory:'
command cd $dir
cd $dir

comment '#put the fortran namelist file in the correct directory'
comment '#to generate some observations for the test:'
command cp make_obs.dat run/pf_parameters.dat
cp make_obs.dat run/pf_parameters.dat

comment '#move to the run directory:'
command cd run
cd run

comment 'Run a single ensemble member to generate observations: '
command /usr/bin/mpirun --output-filename truth -np 1 $mdlex : -np 1 empire
/usr/bin/mpirun --output-filename truth -np 1 $mdlex : -np 1 empire
echo 'Observations generated successfully, look...'
command 'ls obs*'
ls obs*

comment '#move to the 4denvar directory:'
command cd ../../../src/4dEnVar/
cd ../../../src/4dEnVar/

comment '#build 4dEnVar:'
make clean
command make
make > /dev/null
echo '4DEnVar successfully built'


comment '#copy the 4dEnVar executable to the run directory:'
command cp $empex $dir/run/
cp $empex $dir/run/

comment '#move back to the example directory:'
command cd $dir
cd $dir

comment '#copy the 4dEnVar controlling Fortran namelist file'
comment '#to the run directory:'
command cp run_4denvar_cg.nml run/vardata.nml
cp run_4denvar_cg.nml run/vardata.nml

comment '#move to the run directory:'
command cd run
cd run

comment '#run an ensemble of 10 members:'
command /usr/bin/mpirun --output-filename cg -np 24 $mdlex : -np 1 $empex
/usr/bin/mpirun --output-filename cg -np 24 $mdlex : -np 1 $empex

comment '#move back to the original directory:'
command cd $dir
cd $dir

comment '#copy the 4dEnVar controlling Fortran namelist file'
comment '#to the run directory, for lbfgs:'
command cp run_4denvar_lbgfs.nml run/vardata.nml
cp run_4denvar_lbfgs.nml run/vardata.nml

comment '#move to the run directory:'
command cd run
cd run

comment '#run an ensemble of 10 members:'
command /usr/bin/mpirun --output-filename lbfgs -np 24 $mdlex : -np 1 $empex
/usr/bin/mpirun --output-filename lbfgs -np 24 $mdlex : -np 1 $empex

comment '#move back to the original directory:'
command cd $dir
cd $dir



comment '#put the original model specific file back in place'
command cp model_specific_backup.f90 ../../model_specific.f90
cp model_specific_backup.f90 ../../model_specific.f90
