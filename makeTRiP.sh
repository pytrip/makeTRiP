#!/bin/bash
#
# Preparing runs of SHIELD-HIT12A for a range for energies
# for the determination of DDD and SPC files, that is 
# depth-dose distributions and particle type resolved
# spectral files. 
#
# In order to convert the obtained results in the .ddd 
# and .spc formats used by TRiP98 the script postprocess.sh
# has to be run after all output files have been created.
# 
#
# USAGE:   ./makeTRiP.sh <directory>
# EXAMPLE: ./makeTRiP.sh C12lin
# 
# AUTHORS: 
# Niels Bassler <bassler@phys.au.dk>
# Armin Luehr
# Ricky Teiwes
# Toke Printz Ringbaek
#
# DATE: 17.03.2015
#
############################################
#

# read external config file.
source config.sh


# CREATE SUBDIRECTORY
if (( $# < 1)); then
    echo "$0 subdir"
    exit
else
    SUBDIR=$1
    mkdir -p $SUBDIR
fi


BDIR=$PWD

# RANDOM NUMBERS:
RDMNUM=46234201    # DO NOT CHANGE THIS NUMBER
ADDRDM=2327        # THIS NUMBER MAY BE CHANGED

#
for ENG in $(seq $emin $estep $emax)
do
# PREPARING DIRECTORIES
    cd $BDIR
    CDIR=$SUBDIR/`printf "%03d" $ENG`
    EPAD=`printf "%3d" $ENG`
    mkdir -p $CDIR
    echo $CDIR
# COPY REQUIRED INPUT FILES 
    cp detect.dat mat.dat beam.dat geo.dat $CDIR #SHIELD-HIT input files
    cp Air.txt Water.txt Lucite.txt $CDIR #External stopping power tables
    cp 1drifi3.dat 2drifi6.dat $CDIR #External beam modifier data files (for RiFis) 
    cp rtshield.sh $CDIR
    cd $CDIR
# MODIFY INPUT FILES
    sed -i "s/XXX/$EPAD/g" detect.dat
    sed -i "s/XXX.0/$EPAD.0/g" beam.dat
    RDMNUM=$(( $RDMNUM + $ADDRDM ))
    sed -i "s/86732301/$RDMNUM/g" beam.dat
    # sed -i "s/nucdummyfile/            /g" beam.dat 

# START JOB ON CLUSTER
#
# FOR CONDOR
    #nohup rcshield.py -f Polymethyl.txt &

# FOR TORQUE
    qsub -V -d . rtshield.sh              # be mean
#   qsub -V -d . -l nice=19 rtshield.sh   # be nice
done

