The HadCM3/Particle Filter software has two parts, one for the
model and one for the filter. These have to be build separately
but are run in the same script.

Currently it has been set up for run 7 instances of HadCM3, each running on 32
processors.
Each model sends the U and V winds to the PF, then
waits for the perturbed values to be returned and then runs the next timestep.

I've set up various communicators in MPI to allow the transfer of
data.

The PF code hasn't been tested, all that currently happens is that the
U and V winds are returned unaffected.

The PF and model are now completely separate, the PF acts as a "coupler"
for the HadCM3 ensemble. Each model instance runs on it's own node
(32 PEs) and the PF runs in serial on a single node. This means all
of the node memory (96GB?) is available to the PF. In the original code
the PF was serial and I haven't attempted to parallelise it.


As you need the PGI compiler to build the Particle Code, you need to
set this up with

module swap PrgEnv-cray PrgEnv-pgi
module load xt-libnagfl


HadCM3.

To build be model:

Copy xhxrc, user simon in the umui on puma and change the user details
and then submit. There is an extra mod file
/home/n02/n02/simon/pf/test/pf_coupling.mod
where all of the changes required for the filter are set up. To make
futher changes copy this mod and update the umui to use your new version.


PF.

The PF code is located in /home/n02/n02/simon/pf/test. The structure is
different from the orginal version. Use the Makefile to build.
Remember to do

module swap PrgEnv-cray PrgEnv-pgi
module load xt-libnagfl

before running make.
The main control code is in pf_couple.f90 and the MPI routines are in
comms.f90.

Running.

To run the code is a little complicated as the standard UM environment
expects a single model instance controlled by a single script. Running
an ensemble of models is very complex, so I've simplified it for the
development stage.

To set things up:

copy the directory /work/n02/n02/simon/pftmp to your work disk.
copy /work/n02/n02/simon/xhxrb/pf_run.sub to the directory where your
model was built.

Edit pf_run.sub and change the environment variables
TMPDIR to yout pftmp directory
DATAW to your model built/work directory
RUNID to the runid o your model.

Copy pf_couple from your build directory into your new $DATAW

You will have to change the "PBS -A" arguement to your Hector account
group.


Comments:

The initial MPI set up of the Particle Filter is done in  initialise_mpi.
This creates a number of communicators. The one you need to add extra
MPI calls is COUPLE_MPI_COMMUNICATOR. Do not use MPI_COMM_WORLD.

The HadCM3 atmosphere prognotic variables are on different grids. PSTAR
is a single level field, while U, V, QT and THETA_L are on 17 levels.
Also U and V have one fewer grid box in the NS direction (Arakawa B grid).
I don't know if this will affect the analysis part of the PF.

Currently only U and V winds are shared.

Currently I've set the arrays up to use the full atmosphere grid. The final
row of U and V are filled with zeros.

I've used PsiGrand as the main array to hold the data from HadCM3. It has
three dimentions:

nxn*nyn,fields,ensemble

where nxn and nyn are the horizontal dimentions, fields is the total number
of verticle fields and ensemble is the number of model instances. Currently
fields is set to be 34 (2*17) as only the U and V winds are shared between
the PF and model.

The form of PsiGrand isn't fixed, it can be easily changed, if required.

The analsys stage of the code hasn't been tested, and I'll be very
surprised if it works with the new model. It is compiled but not called.

The nudging part of the PF exists as skeleton code, but does nothing
at the moment.

The models are currently configured to run for a day.

The submission system is a little inellegant, espically the model
configuration which currently can't be changed (things like run length).
I will set up a generic system in due course.

There's likely to be an issue with the diagnostic output from the model.
Currently it's possible that the inidivual models will overwrite each
other's output data as they will attempt to write to files of the same 
filename. This is being investigated.

The output is currently riddled with debug print statements, feel
free to remove them.

We will need a system to develope code in parallel. I'm minded to give
git a try, as this is the preferred code management system currently.

















