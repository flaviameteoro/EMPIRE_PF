#!/bin/bash
set -o verbose

#make the lorenz96 code
make lorenz96

#backup the model_specific.f90 file
cp model_specific.f90 model_specific.f90_original

#make all the changes to model_specific.f90 that we have listed above
#happily they are in the file examples/lorenz96/model_specific_tutorial_l96.f90
cp examples/lorenz96/model_specific_tutorial_l96.f90 model_specific.f90

#build the codes
make

#check to see if the empire codes built properly
ls -l bin/empire

#make a run directory
mkdir -p rundirectory

#move to the run directory
cd rundirectory

#get the empire.nml file
cp ../examples/lorenz96/tutorial1.nml empire.nml

#look at the empire.nml file
cat empire.nml

#pause to look at this file
sleep 10

#generate the l96.nml file
echo -e "&l96\nN=40,\ntotal_timesteps=4,\nF=8.0d0,\ndt=1.0d-2\n/\n" > l96.nml

#look at the l96.nml file
cat l96.nml

#pause to look at this file
sleep 5

#generate the observations
mpirun --output-filename truth -np 1 ../bin/lorenz96 : -np 1 ../bin/empire

#look for the observation files
ls obs*

#modify the empire.nml file to run a stochastic ensemble
sed -i "s/filter.*/filter='SE',/g" empire.nml
sed -i "/gen_data/d" empire.nml

#look at the empire.nml file
cat empire.nml

#pause to look at this file
sleep 10

#now run the stochastic ensemble 
mpirun --output-filename stoch -np 32 ../bin/lorenz96 : -np 4 ../bin/empire

#modify the empire.nml file to run the LETKF
sed -i "s/filter.*/filter='LE',/g" empire.nml

#look at the empire.nml file
cat empire.nml

#pause to look at this file
sleep 5

#now run the LETKF 
mpirun --output-filename assim -np 32 ../bin/lorenz96 : -np 4 ../bin/empire

#plot the output
../examples/lorenz96/tutorial_lorenz96_plot.py
