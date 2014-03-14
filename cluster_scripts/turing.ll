#!/bin/bash

# submit with llsubmit
# only for use with Turing.  Does not yet work with restart.

#=========== Global directives ===========
#@ job_name = unit_testing
#@ shell    = /bin/bash

#============= pre-parallel step =============
#@ step_name = sequential_preprocessing
#@ job_type  = serial
#@ cpu_limit = 0:02:00
#@ output = $(job_name).$(jobid).pre
#@ error = $(output)
#@ queue

#============= Parallel step =============
#@ step_name  = parallel_step
#@ dependency = (sequential_preprocessing == 0)
# (submit only if previous step completed without error)
#@ job_type = BLUEGENE
#@ wall_clock_limit = 00:10:00
#@ bg_size = 64
#@ output = $(job_name).$(jobid)
#@ error = $(output)
#@ queue

#============= post-parallel step =============
#@ step_name  = sequential_postprocessing
#@ dependency = (parallel_step >= 0)
# (submit even if previous step completed with an error)
#@ job_type   = serial
#@ cpu_limit  = 0:02:00
#@ output = $(job_name).$(jobid).post
#@ error = $(output)
#@ queue

PARAMSFILE=Testing_Sphere.ini
PROGFILE=flusi
# INPUT FILES FOR THE SIMULATION
files_input=($PARAMSFILE $PROGFILE)
# when retaking a backup, the *.t files as well as runtime_backupX.h5
# must be copied


case $LOADL_STEP_NAME in

  #======= Sequential preprocessing ========
    sequential_preprocessing )
	set -ex
	
	echo "start dir:" $LOADL_STEP_INITDIR

	if [ ! -f  $LOADL_STEP_INITDIR/$PROGFILE ]; then
	    echo  $LOADL_STEP_INITDIR/$PROGFILE "not found!"
	    exit
	fi

	if [ ! -f $PARAMSFILE ]; then
	    echo  $LOADL_STEP_INITDIR/$PARAMSFILE "not found!"
	    exit
	fi

	# copy all input files to the work directory
        for file in "${files_input[@]}"
        do
          cp $LOADL_STEP_INITDIR/$file $tmpdir  
        done

	echo "working dir: " $tmpdir
	cd $tmpdir
	
	echo "files in starting dir:"
	ls -l
	;;

  #============= Parallel step =============
    parallel_step )
	set -x
	cd $tmpdir

	pwd

	ls -l

 	runjob --np 1024 --ranks-per-node 16 : ./$PROGFILE $PARAMSFILE
	;;

  #======= Sequential postprocessing =======
    sequential_postprocessing )
	set -x

	mv $tmpdir/* $LOADL_STEP_INITDIR/

	;;
esac
