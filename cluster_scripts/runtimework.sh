#!/bin/bash

# copy the runtime backup files files from working directory on turing
# to pwd

# usage: ./runtimework.sh <run directory on turing>

# example:
# ./runtimework.sh /linkhome/rech/sch/rsch511/flusi/smc/bc256/qay1/70.6/workdir/
# note that the work dir must have a trailing slash.

if [ "$1" == "" ]; then
    srcdir=$(pwd | sed 's/^.*lidris//')
    srcdir=/linkhome/rech/sch/rsch511/flusi${srcdir}/
else
    srdcir=$1
fi

echo $srcdir

ionice -c 3 nice -n 19 rsync -avpun  --include 'runtime*.h5' --exclude '*' turing:${srcdir} ./

echo "pull files? Y/N"
read text

if [ "$text" == "Y" ] ; then
    ionice -c 3 nice -n 19 rsync -avpu  --include 'runtime*.h5' --exclude '*' turing:${srcdir} ./
fi
