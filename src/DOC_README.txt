/** \mainpage EMPIRE Data Assimilation Documentation
\author Philip A. Browne
\date Time-stamp: <2014-09-23 12:07:45 pbrowne>

\section Download Downloading

These codes are hosted on www.bitbucket.org and can be attained with the following commands:
\code{.sh}
git clone git@bitbucket.org:pbrowne/empire-data-assimilation.git
\endcode

or

\code{.sh}
wget https://bitbucket.org/pbrowne/empire-data-assimilation/get/aa31fdfc3912.zip && gunzip aa31fdfc3912.zip
\endcode

\section Compiling Compiling

\subsection Compiling_code Compilation of the source code

The Makefile must be editted for the specific compiler setup. In the main directory you will find the file \c Makefile.

Edit the variables as follows:
- \c FC The fortran compiler
- \c FCOPTS The options for the fortran compiler
- \c LIB_LIST The libraries to be called. Note this must include BLAS

To compile the source code, simply then type the command
\code{.sh}
make
\endcode

If successful, the following executables are created in the bin/ folder:
- \link empire \endlink
- \link alltests \endlink
- \link test_h \endlink
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

- \link h \endlink This is the observation operator
 
- \link ht \endlink This is the transpose of the observation operator

- \link r \endlink This is the observation error covariance matrix \f$R\f$

- \link rhalf \endlink This is the square root of the observation error covariance matrix \f$R^{\frac{1}{2}}\f$
 
- \link solve_r \endlink This is a linear solve with the observation error covariance matrix, i.e. given \f$b\f$, find \f$x\f$ such that \f$Rx=b\f$ or indeed, \f$x = R^{-1}b\f$

- \link q \endlink This is the model error covariance matrix \f$Q\f$

- \link qhalf \endlink This is the square root model error covariance matrix \f$Q^{\frac{1}{2}}\f$
 
- \link solve_hqht_plus_r \endlink This is a linear solve with the matrix \f$(HQH^T+R)\f$

Not all of these subroutines will be required for each filtering method you wish to use, so it may be advantageous to only implement the necessary ones.

\section Testing Testing 

You can test your user supplied routines by running the test codes found in the folder bin/.

These are by no means full-proof ways of ensuring that you have implemented things correctly, but should at least check what you have done for logical consistency.

For example, they will test if \f$ H H^Tx = x\f$, and if \f$ Q^{\frac{1}{2}}Q^{\frac{1}{2}}x = Qx\f$ for various different vectors \f$x\f$.

\section Linking Linking to your model using EMPIRE
\todo Write a bit about how to put the MPI commands into the model, or point to where that is. <a href="http://www.met.reading.ac.uk/~darc/empire">www.met.reading.ac.uk/~darc/empire</a>

\section Running Running
For example, to run \b N_MDL copies of the model with \b N_DA copies of empire, then the following are possible:
\code{.sh}mpirun -np N_MDL model_executable : -np N_DA empire\endcode

\code{.sh}aprun -n N_MDL -N N_MDL model_executable : -n N_DA -N N_DA empire\endcode

\todo I have to talk about how pf_parameters.dat works.






\section Bugs Bug Reports and Functionality Requests
While the code is not too large, you may email me the issue or request <a href="mailto:p.browne@reading.ac.uk">here</a>.


However there is a webpage set up for this:

<a href="https://bitbucket.org/pbrowne/empire-data-assimilation/issues">https://bitbucket.org/pbrowne/empire-data-assimilation/issues</a>


*/
