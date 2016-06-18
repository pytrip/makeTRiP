#!/bin/sh
#
# how to use this
# change to directory you want to run
# $ qsub -V -t 0-9 -d . rtshield.sh
#
#PBS -N SHIELDHIT_JOB
#PBS -l walltime=12:00:00

#print the time and date
#date

#echo Task number $PBS_TASKNUM
#echo Node number $PBS_NODENUM
#echo Job id      $PBS_O_JOBID
#echo Array id    $PBS_ARRAYID
#echo $PWD
#echo $PATH

shieldhit $PBS_ARRAYID

