#!/bin/bash
#
# Performing postprocessing of data files produced by SHIELD-HIT12A 
# using the script makeTRiP.sh. DDD and SPC files, that is 
# depth-dose distributions and particle type resolved spectral 
# files, are produced in the .ddd and .spc formats used by TRiP98.
#
#
# USAGE:   ./postprocess.sh <subdir>
# EXAMPLE: ./postprocess.sh Carbon_RiFi2D6
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
#
CLEAND=N

# read external config file.
source config.sh

############################################


if (( $# < 1)); then
    echo "$0 subdir"
    exit
else
    SUBDIR=${1%/}  # Remove possible slash.
fi


BDIR=$PWD
cd $SUBDIR
mkdir -p SPC



####
# convert and move output files from simulation
#

for ENG in $(seq $emin $estep $emax)
do

    cd $BDIR
    CDIR=$SUBDIR/`printf "%03d" $ENG`
    EPAD=`printf "%3d" $ENG`
    E3=`printf "%03d" $ENG`
    cd $CDIR
    echo $EPAD 

# convert 2D DDD raw data file to ascii
#    if [ -e $ENG"CrD" ]; then
      bdo2txt "$ENG"CrD0000.bdo #changed to bdo2txt and file name accordingly /TPR 
      mv "$ENG"CrD0000.txt ../DDD_Cr_"$E3".dat #changed file name to new output format /TPR
#    fi
#TODO: maybe change names in detect.dat? Cr is for carbon ions so here it should be Hr? But required further changes in ddd.py /TPR

# convert 1D DDD raw data file to ascii
#    if [ -e $ENG"CrD1" ]; then
      bdo2txt "$ENG"CrD10000.bdo #changed to bdo2txt and file name accordingly /TPR 
      mv "$ENG"CrD10000.txt ../DDD_Cr1_"$E3".dat #changed file name to new output format /TPR
#    fi

# Copy SPC file respecting the TRiP naming convention
# NBassler: COPY don't move, since move will destroy the data source.    
    cp "$ENG"Cr.spc0000.bdo        ../SPC/${ion}.H2O.MeV"$E3"00.spc #changed file name to new output format /TPR

done

####
# run ddd.py script in order to produce .ddd files in directory DDD
#

cd $BDIR
cd $SUBDIR
python ../ddd.py 

#### CLEAN SUBDIRECTORIES
# 

if [ "$CLEAND" == "Y" ]; then
  echo "Cleaning subdirectories from temporary files"
  for ENG in $(seq $emin $estep $emax)
  do

    cd $BDIR
    CDIR=$SUBDIR/`printf "%03d" $ENG`
    EPAD=`printf "%3d" $ENG`
    cd $CDIR

    rm for0*00* 
    cd ..
  done
fi

####
# build .tar.gz files for SPC and DDD data
#
echo Build SPC_${SUBDIR}.tar.gz
tar -czf SPC_${SUBDIR}.tar.gz SPC/
echo Build DDD_${SUBDIR}.tar.gz #added for DDD similar line as for SPC /TPR
tar -czf DDD_${SUBDIR}.tar.gz DDD/
