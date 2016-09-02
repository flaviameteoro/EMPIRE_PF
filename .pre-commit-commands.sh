#!/bin/bash
set -o errexit
#COMMANDS THAT ARE RUN BEFORE A GIT COMMIT IS MADE


#mv to the home folder:
cd $(git rev-parse --show-toplevel)

#run the compilation tests
#./tests/test_compile.sh
