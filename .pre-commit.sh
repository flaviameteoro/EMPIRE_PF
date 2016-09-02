#!/bin/bash

#stash the things which arent commited
git stash -q --keep-index
#run the commands that we want
./.pre-commit-commands.sh
RESULT=$?
#put the non-commited changes back
git stash pop -q
# if the commands failed, then exit and do not commit
[ $RESULT -ne 0 ] && exit 1
exit 0
