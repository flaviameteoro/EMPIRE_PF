#!/bin/bash
set -o nounset
set -o errexit

script=$(basename "$0")
function is_integer() { [ "$1" -eq "$1" ] > /dev/null 2>&1 ; return $? ; }
function usage() {
    cat <<EOF
$script usage :::

    $script takes a single argument that corresponds to the
    number of ensemble members to be used.

    Examples:
    1) Running a single ensemble member
          $script 1

    2) Running an ensemble of 20 members
          $script 20 

    PAB 9 June 2015
EOF
}
function return_error() { usage ; exit -1 ; }
function echor() { echo $* ; $* ; }
function eecho() { echo $* >&2 ; return_error ; }
function jecho() { echo $* >> $js ; }

if [ "$#" == 0 ] || [ "$#" -gt 1 ] ; then
    eecho "ERROR: Incorrect number of arguments to $script"
fi

if ! is_integer $1; then
    eecho "ERROR: Argument to $script not an integer"
fi

empire_ex=../bin/empire_linear_identity
model_ex=../bin/linear_empire_vader
namelist=empire.nml
mpirun=mpirun
js=jobscript.sh
wd=$PWD


#test for namelist
if [ ! -f $namelist ]; then
    eecho "Error:: no namelist $namelist found"
fi

#define output names
if grep --quiet "gen_data=.true." $namelist || grep --quiet "gen_data=T" $namelist; then
    runname=truth
else
    filter=$(grep filter $namelist | cut -f2 -d= | cut -f2 -d"'")
    echo selected filter = $filter
    runname=$filter
fi


echo Running an ensemble of $1 members
echo Output will be prefixed with the name $runname


cat <<EOF > $js
#!/bin/bash
$mpirun --output-filename $runname -np $1 $model_ex : -np 1 $empire_ex
EOF
chmod 744 $js


echo "The following jobscript, $js, will be submitted for batch processing:"
echo -e "\n \n \n \n"
cat $js
echo -e "\n \n \n \n"

#replace the line below with a suitable submission line on your own system
#my testing system has no queue, so I just run the jobscript directly:
echor $js

# on the UoR Metcluster, and example would be the following:

#let np=$1+1
#np=$(($np<8?$np:8))
#echor qsub -l hpitt,limited -R y -pe openmpi $np $js

