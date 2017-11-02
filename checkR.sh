#!/bin/bash
set -e

exitstatus=0

Rscript -e "lintr::lint(PRSice.R)"
outputbytes=`Rscript -e "lintr::lint(PRSice.R)" | grep ^ | wc -c`
if [ $outputbytes -gt 0 ]
then
  exitstatus=1
fi

exit $exitstatus
