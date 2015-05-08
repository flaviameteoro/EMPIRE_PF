/** \mainpage EMPIRE Data Assimilation Documentation
\author Philip A. Browne <a href="mailto:p.browne@reading.ac.uk">p.browne@reading.ac.uk</a>
\date Time-stamp: <2015-05-08 10:57:22 pbrowne>

\section Methods EMPIRE Methods
For a list of methods implemented in EMPIRE, please click here: \link methods \endlink

\section Download Downloading

The current version is an *academic version*, and the user interface may be subject to change.

These codes are hosted on www.bitbucket.org and can be obtained with the following commands:
\code{.sh}
git clone https://www.bitbucket.org/pbrowne/empire-data-assimilation.git
\endcode
To upgrade to the latest versions of the codes, use the following command:
\code{.sh}
git pull https://www.bitbucket.org/pbrowne/empire-data-assimilation.git
\endcode

\copyright These codes are distributed under the GNU GPL v3 License. See LICENSE.txt.

\section Compiling Compiling

\subsection Compiling_code Compilation of the source code

The Makefile must be edited for the specific compiler setup. In the main directory you will find the file \c Makefile.in.

Edit the variables as follows:
- \c FC The fortran compiler

This has been tested with gfortran 4.8.2, crayftn 8.2.6 and ifort 14.0.1.106

- \c FCOPTS The options for the fortran compiler
- \c LIB_LIST The libraries to be called. Note this must include BLAS and LAPACK

- \c MODFLAG The flag to specify where module files should be placed by the fortran complier. Examples are
  - \c gfortran: -J
  - \c ifort: -module
  - \c crayftn: -em -J
  - \c pgfortran: -module

To compile the source code, simply then type the command
\code{.sh}
make
\endcode

If successful, the following executables are created in the bin/ folder:
- \link empire \endlink
- \link alltests \endlink
<!-- - \link test_h \endlink -->
- \link test_hqhtr \endlink
- \link test_q \endlink
- \link test_r \endlink

To remove the object and executable files if compilation fails for some reason, run the following:
\code{.sh}make clean\endcode

\subsection Compiling_docs Compilation of the documentation
Documentation of the code is automatically generated using Doxygen, dot and pdflatex.

All of these packages must be installed for the following to work.

\code{.sh}make docs\endcode

This will make an html webpage for the code, the mainpage for which is located in doc/html/index.html.

A latex version of the documentation will be built to the file doc/latex/refman.pdf.

To simply make the html version of the documentation (if pdflatex is not available) then use the command \code{.sh}make doc_html\endcode


\section Custom Customising for specific models

<em>This is where the science and all the effort should happen!!</em>

The file model_specific.f90 should be editted for the specific model which you wish to use. This contains a number of subroutines which need to be adapted for the model and the observation network. We list these subsequently.

- \link configure_model \endlink This is called early in the code and can be used to read in any data from files before subsequently using them in the below operations.

- \link reconfigure_model \endlink This is called after each observation timestep. If the observation dimension changes it should be updated here, along with the number of model timesteps until the next observation

- \link h \endlink This is the observation operator
 
- \link ht \endlink This is the transpose of the observation operator

- \link r \endlink This is the observation error covariance matrix \f$R\f$

- \link rhalf \endlink This is the square root of the observation error covariance matrix \f$R^{\frac{1}{2}}\f$
 
- \link solve_r \endlink This is a linear solve with the observation error covariance matrix, i.e. given \f$b\f$, find \f$x\f$ such that \f$Rx=b\f$ or indeed, \f$x = R^{-1}b\f$

- \link solve_rhalf \endlink This is a linear solve with the square root of the observation error covariance matrix, i.e. given \f$b\f$, find \f$x\f$ such that \f$R^{\frac{1}{2}}x=b\f$ or indeed, \f$x = R^{-\frac{1}{2}}b\f$

- \link q \endlink This is the model error covariance matrix \f$Q\f$

- \link qhalf \endlink This is the square root model error covariance matrix \f$Q^{\frac{1}{2}}\f$
 
- \link solve_hqht_plus_r \endlink This is a linear solve with the matrix \f$(HQH^T+R)\f$

- \link dist_st_ob \endlink This specifies the distance between a an element of the state vector and an element of the observation vector

- \link bhalf \endlink This is the square root of the background error covariance matrix \f$B^{\frac{1}{2}}\f$

- \link get_observation_data \endlink This subroutine must return the observation data at, or subsequently to, the given timestep. This routine only needs to be editted if you wish to use your own observations. It is set up to work automatically with pseudo-observations for running twin experiments.

Not all of these subroutines will be required for each filtering method you wish to use, so it may be advantageous to only implement the necessary ones.

\section Testing Testing 

You can test your user supplied routines by running the test codes found in the folder bin/.

These are by no means full-proof ways of ensuring that you have implemented things correctly, but should at least check what you have done for logical consistency.

For example, they will test if \f$ R^{-1} Ry = y\f$, and if \f$ Q^{\frac{1}{2}}Q^{\frac{1}{2}}x = Qx\f$ for various different vectors \f$x, y\f$.

\section Linking Linking to your model using EMPIRE
Full instructions on how to put the EMPIRE MPI commands into a new model can be found at <a href="http://www.met.reading.ac.uk/~darc/empire">www.met.reading.ac.uk/~darc/empire</a>.

\section Running Running
For example, to run \b N_MDL copies of the model with \b N_DA copies of empire, then the following are possible:
\code{.sh}mpirun -np N_MDL model_executable : -np N_DA empire\endcode

\code{.sh}aprun -n N_MDL -N N_MDL model_executable : -n N_DA -N N_DA empire\endcode

The empire executable is controlled by the namelist data file \link pf_control::parse_pf_parameters pf_parameters.dat\endlink. As such, this file should be put in the directory where empire is executed.

\section Examples Examples
In the directory \c examples there is currently one example of how to use EMPIRE, specifically with the Lorenz 1996 model. In the directory you will find an example model_specific.f90 file setup for that model, along with a file \c instructions.txt which will lead you step by step through how to run a twin experiment.




\section Bugs Bug Reports and Functionality Requests
While the code is not too large, you may email me the issue or request <a href="mailto:p.browne@reading.ac.uk">here</a>.


However there is a webpage set up for this:

<a href="https://bitbucket.org/pbrowne/empire-data-assimilation/issues">https://bitbucket.org/pbrowne/empire-data-assimilation/issues</a>


*/

/*! \page methods Assimilation Methods
\section methods_filters Filters
The filters implemented in EMPIRE can be divided into two categories, particle filters and Ensemble Kalman filters

\subsection methods_pfs Particle filters
\subsubsection SIRfilter SIR Filter (Sequential Importance Resampling) 
See file @ref sir_filter \n
<a href="http://dx.doi.org/10.1049/ip-f-2.1993.0015">Gordon, Salmond and Smith (1993)</a>.\n
Model specific operations required: \n
  - \link qhalf \endlink
  - \link h \endlink
  - \link solve_r \endlink \n
The SIR filter has no parameters to be chosen. \n
To select the SIR filter, in \link pf_control::parse_pf_parameters pf_parameters.dat \endlink set the following variables:
   - \link pf_control::pf_control_type::filter filter \endlink = 'SI'
\subsubsection EWPF Equivalent Weights Particle Filter 
See files @ref proposal_filter @ref equivalent_weights_filter\n
<a href="http://doi.wiley.com/10.1002/qj.699">Van Leeuwen (2010)</a>.\n
Model specific operations required:
 - \link qhalf \endlink
 - \link q \endlink
 - \link h \endlink
 - \link ht \endlink
 - \link solve_r \endlink
 - \link solve_hqht_plus_r \endlink
 - \link rhalf \endlink \n
The Equivalent Weights particle filter has a number of free parameters to be chosen. \n
   - \link pf_control::pf_control_type::nudgefac nudgefac\endlink
   - \link pf_control::pf_control_type::nfac nfac \endlink
   - \link pf_control::pf_control_type::ufac ufac \endlink
   - \link pf_control::pf_control_type::keep keep \endlink \n
To select the EWPF, in \link pf_control::parse_pf_parameters pf_parameters.dat \endlink set the following variables:
   - \link pf_control::pf_control_type::filter filter \endlink = 'EW'
\subsection methods_enkfs Ensemble Kalman filters

\subsubsection LETKF LETKF (The Localised Ensemble Transform Kalman Filter)
See file  @ref letkf_analysis \n
<a href="http://dx.doi.org/10.1016/j.physd.2006.11.008">Hunt, Kostelich and Szunyogh (2007)</a>. \n
Model specific operations required: \n
 - \link h \endlink
 - \link solve_rhalf \endlink
 - \link dist_st_ob \endlink \n
The LETKF has a number of free parameters to be chosen. \n
   - \link pf_control::pf_control_type::rho rho \endlink
   - \link pf_control::pf_control_type::len len \endlink \n 
To select the LETKF, in \link pf_control::parse_pf_parameters pf_parameters.dat \endlink set the following variables:
   - \link pf_control::pf_control_type::filter filter \endlink = 'LE' or 'LD' with LE including model error but LD being deterministic


\section methods_smoothers Smoothers

Coming at some point in the future: LETKS (Please contact us if you want us to develop this sooner rather than later)


\section methods_var Variational Methods

\subsection fourdEnVar 4dEnVar

See files @ref fourdenvar @ref fourdenvar_fcn \n
<a href="http://dx.doi.org/10.1175/2008MWR2312.1">Liu, Xian and Wang (2008)</a>.\n

Model specific operations required:
 - \link h \endlink
 - \link solve_r \endlink

Currently there is the basic functionality to do 4dEnVar so long as EMPIRE-VADER is used for reverse communication. This is work in progress.

\todo Add some stuff about how to use this.

*/


/*! \page features Other EMPIRE features
\section gendata Generating artificial observations

EMPIRE can generate artificial observations easily and quickly.

Model specific operations required: \n
 - \link h \endlink
 - \link rhalf \endlink
 - \link qhalf \endlink

In \link pf_control::parse_pf_parameters pf_parameters.dat \endlink set the following variables:
- \link pf_control::pf_control_type::gen_data gen_data \endlink = .true.

The system then should be run with a single ensemble member and a single EMPIRE process, i.e.
\code{.sh}
mpirun -np 1 model : -np 1 empire
\endcode

\section Observations Observations
To use real observations (i.e. those not generated automatically in twin experiment mode) the user must change the subroutine \link get_observation_data get_observation_data\endlink in model_specific.f90.

When called, \link get_observation_data get_observation_data\endlink must return the vector of observations \f$y\f$ that corresponds to the observation on, subsequently to, the current timestep.

\section detens Running a deterministic ensemble

EMPIRE can simply integrate forward in time an ensemble of models.

In \link pf_control::parse_pf_parameters pf_parameters.dat \endlink set the following variables:
- \link pf_control::pf_control_type::filter filter \endlink = 'DE'

\section stochens Running a stochastic ensemble

EMPIRE can integrate forward in time an ensemble of models whilst adding stochastic forcing.

Model specific operations required: \n
 - \link qhalf \endlink

In \link pf_control::parse_pf_parameters pf_parameters.dat \endlink set the following variables:
- \link pf_control::pf_control_type::filter filter \endlink = 'SE'

\section rankhistograms Outputting rank histograms

This is controlled by \link pf_control::pf_control_type::use_talagrand use_talagrand \endlink in \link pf_control::parse_pf_parameters pf_parameters.dat \endlink  and for more information see \link histogram_data::load_histogram_data load_histogram_data \endlink.

\section trajectories Outputting trajectories of model variables

This is controlled by \link pf_control::pf_control_type::use_traj use_traj \endlink in \link pf_control::parse_pf_parameters pf_parameters.dat \endlink  and for more information see \link traj_data::setup_traj setup_traj \endlink.
*/

/*! \page citing How to Cite EMPIRE

# EMPIRE itself

For all applications that use these codes, please cite the following paper:

PA Browne, S Wilson (2015)

A simple method for integrating a complex model into an ensemble data assimilation system using MPI

<http://dx.doi.org/10.1016/j.envsoft.2015.02.003>

# Use of different methods within EMPIRE

## Equivalent weights particle filter
Van Leeuwen (2010)

Nonlinear data assimilation in geosciences: an extremely efficient particle filter

<http://doi.wiley.com/10.1002/qj.699>

## Sequential importance resampling
Gordon, Salmond and Smith (1993)

Novel approach to nonlinear/non-Gaussian Bayesian state estimation

<http://dx.doi.org/10.1049/ip-f-2.1993.0015>

## Localised Ensemble Transform Kalman Filter
Hunt, Kostelich and Szunyogh (2007)

Efficient data assimilation for spatiotemporal chaos: A local ensemble transform Kalman filter

<http://dx.doi.org/10.1016/j.physd.2006.11.008>

## 4DEnVar
Liu, Xian and Wang (2008)

An Ensemble-Based Four-Dimensional Variational Data Assimilation Scheme. Part I: Technical Formulation and Preliminary Test

<http://dx.doi.org/10.1175/2008MWR2312.1>

# Use of different external codes within EMPIRE

## CG+
Gilbert and Nocedal (1992)

Global Convergence Properties of Conjugate Gradient Methods for Optimization

<http://dx.doi.org/10.1137/0802003>

Software available here: <http://users.iems.northwestern.edu/~nocedal/CG+.html>


## L-BFGS-B 

Byrd, Lu and Nocedal (1995)

A Limited Memory Algorithm for Bound Constrained Optimization

<http://dx.doi.org/10.1137/0916069>

and/or

Zhu, Byrd and Nocedal (1997)

L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization

<http://dx.doi.org/10.1145/279232.279236>

Software available here: <http://users.iems.northwestern.edu/~nocedal/lbfgsb.html>


*/
