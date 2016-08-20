#!/bin/bash
if [ $# -ne 1 ]
then
  echo "... ERROR: specify the commit message."
  exit 1
fi
#
git add *
git status
git commit -a -m "$1"
git remote add github git@github.com:jeokim/camille.git
git push github dev
