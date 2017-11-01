#!/bin/bash
if [ $# -ne 1 ]
then
  echo "... WARNING: specify a commit message."
  MSG="updated"
else
  MSG=$1 
fi
#
git add *
git status
git commit -a -m "$MSG"
git remote add github git@github.com:jeokim/camille.git
git push github LEE_linearNozzle
