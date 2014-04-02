#!/bin/bash

# copy xmf and h5 files from working directory on turing to pwd

srcdir=$(pwd | sed 's/^.*lidris//')
srcdir=/linkhome/rech/sch/rsch511/flusi${srcdir}/workdir/

echo $srcdir

ionice -c 3 nice -n 19 rsync -avpu --remove-source-files --exclude='runtime*' --include='*/' --include='*.xmf' --include='*.h5' --exclude='*'  turing:${srcdir} ./

#ionice -c 3 nice -n 19 rsync -avpun turing:${srcdir} ./
