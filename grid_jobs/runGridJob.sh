#!/usr/bin/env bash
# runGridJob.sh

# Expects a tar file that unpacks to lib/ and fcl/ directories

# Get things started!
date
hostname
echo $CLUSTER
echo $PROCESS
echo $JOBSUBJOBID
pwd

# Determine location of unpacked libraries
D=${CONDOR_DIR_INPUT}/irmaData
echo "D=$D"
ls $D

# Set up studio environment
source /cvmfs/gm2.opensciencegrid.org/prod/g-2/setup
export PRODUCTS=/cvmfs/gm2.opensciencegrid.org/specials/sl7/prod/g-2:/cvmfs/gm2.opensciencegrid.org/specials/sl7/prod/external:$PRODUCTS 
setup gm2 v9_52_00 -q prof
setup cmake v3_17_1
setup hep_hpc v0_13_01 -q e15
setup ifdh_art v2_05_05 -q e15:prof:s65
export LD_LIBRARY_PATH=$D/lib:$LD_LIBRARY_PATH

# Setup IFDH and SAM
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
setup fife_utils
export SAM_EXPERIMENT=gm2
export EXPERIMENT=gm2

# Get the SAM Project URL (jobsub will pass in the SAM_PROJECT_NAME)
SAM_PROJECT_URL=$(ifdh findProject ${SAM_PROJECT_NAME} gm2 2>&1)

# Establish the SAM process
HOST=`/bin/hostname`
SAM_PROCESS_ID=$(ifdh establishProcess ${SAM_PROJECT_URL} gm2 v9_52_00 $HOST lyon irmaData $JOBSUBJOBID 1 xroot 2>&1)
# Format of above ... ifdh establishProcess url app version location user appfamily description filelimit schema

# Run it!
gm2 -c $D/fcl/irmaData.fcl --sam-web-url=${SAM_PROJECT_URL} --sam-process-id=${SAM_PROCESS_ID}

# Copy the output back
OUT=/pnfs/GM2/scratch/users/lyon/irmaData

mv irmaData.h5 irmaData_${CLUSTER}_${PROCESS}.h5
ifdh cp -D irmaData_${CLUSTER}_${PROCESS}.h5 ${OUT}

mv irmaData.log irmaData_${CLUSTER}_${PROCESS}.log
ifdh cp -D irmaData_${CLUSTER}_${PROCESS}.log ${OUT}
