#!/bin/bash

# copy xmf and h5 files from working directory on turing to pwd

srcdir=$(pwd | sed 's/^.*lidris//')
srcdir=/linkhome/rech/sch/rsch511/flusi${srcdir}/

echo $srcdir

ionice -c 3 nice -n 19 rsync -avpun --remove-source-files --exclude='runtime*' --include='*/' --include='*.xmf' --include='*.h5' --exclude='*'  turing:${srcdir} ./

echo "pull files? Y/N"

read text

if [ "$text" == "Y" ] ; then
    ionice -c 3 nice -n 19 rsync -avpu --remove-source-files --exclude='runtime*' --include='*/' --include='*.xmf' --include='*.h5' --exclude='*'  turing:${srcdir} ./
fi
