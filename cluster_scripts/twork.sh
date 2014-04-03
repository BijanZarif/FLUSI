#!/bin/bash

# copy non-(hdf or xmf) files from the workdir on an active turing run
# to the pwd.

# usage: ./twork.sh <run directory on turing>

# example: 
# ./twork.sh /linkhome/rech/sch/rsch511/flusi/smc/bc256/qay1/70.6/workdir/
# note that the work dir must have a trailing slash.

if [ "$1" == "" ]; then
    srcdir=$(pwd | sed 's/^.*lidris//')
    srcdir=/linkhome/rech/sch/rsch511/flusi${srcdir}/workdir/
else
    srdcir=$1
fi    

echo $srcdir

ionice -c 3 nice -n 19 rsync -avpun --exclude 'core*' --exclude '*.h5' --exclude '*.xmf'  turing:${srcdir} ./

echo "pull files? Y/N"
read text

if [ "$text" == "Y" ] ; then
    ionice -c 3 nice -n 19 rsync -avpu --exclude 'core*' --exclude '*.h5' --exclude '*.xmf'  turing:${srcdir} ./
fi
