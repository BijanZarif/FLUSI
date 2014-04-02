#!/bin/bash

# copy the runtime backup files files from working directory on turing
# to pwd

# usage: ./runtimework.sh <run directory on turing>

# example:
# ./runtimework.sh /linkhome/rech/sch/rsch511/flusi/smc/bc256/qay1/70.6/workdir/
# note that the work dir must have a trailing slash.

if [ "$1" == "" ]; then
    srcdir=$(pwd | sed 's/^.*lidris//')
    srcdir=/linkhome/rech/sch/rsch511/flusi${srcdir}/workdir/
else
    srdcir=$1
fi

echo $srcdir

ionice -c 3 nice -n 19 rsync -avpu --exclude '*' --include 'runtime*.h5' turing:${srcdir} ./
