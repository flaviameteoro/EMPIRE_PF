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
mdlex=linear_empire_vader
empex=empire_4denvar
n=4


comment "#move to the top level directory and build the model:"
command cd ../../
cd ../../
command make linear
make linear
echo '#Model built successfully'

comment '#move back to this example directory:'
command cd $dir
cd $dir

comment '#move the model to a "run folder"'
command mkdir -p run/
mkdir -p run/
command mv ../../bin/$mdlex run/
mv ../../bin/$mdlex run/


comment '#put an appropriate model specific file in the right place:'
command cp -n ../../model_specific.f90 ../../model_specific_backup.f90
cp -n ../../model_specific.f90 ../../model_specific_backup.f90
command cp model_specific_linear.f90 ../../model_specific.f90
cp model_specific_linear.f90 ../../model_specific.f90


comment '#move to the root folder:'
command cd ../../
cd ../../

comment '#Build the full empire code:'
make clean
command make
make > /dev/null 2>&1
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

comment '#move to the root directory:'
command cd ../../../
cd ../../../

comment '#build 4dEnVar:'
command make 4dEnVar
make 4dEnVar > /dev/null
echo '4DEnVar successfully built'


comment '#copy the 4dEnVar executable to the run directory:'
command cp bin/$empex $dir/run/
cp bin/$empex $dir/run/

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

comment "#run an ensemble of $n members:"
command /usr/bin/mpirun --output-filename cg -np $n $mdlex : -np 1 $empex
/usr/bin/mpirun --output-filename cg -np $n $mdlex : -np 1 $empex

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

comment "#run an ensemble of $n members:"
command /usr/bin/mpirun --output-filename lbfgs -np $n $mdlex : -np 1 $empex
/usr/bin/mpirun --output-filename lbfgs -np $n $mdlex : -np 1 $empex

comment '#move back to the original directory:'
command cd $dir
cd $dir



comment '#put the original model specific file back in place'
command mv ../../model_specific_backup.f90 ../../model_specific.f90
mv ../../model_specific_backup.f90 ../../model_specific.f90
