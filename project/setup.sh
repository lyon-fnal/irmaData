#!/usr/bin/env bash
source /cvmfs/gm2.opensciencegrid.org/prod/g-2/setup
export PRODUCTS=/cvmfs/gm2.opensciencegrid.org/specials/sl7/prod/g-2:/cvmfs/gm2.opensciencegrid.org/specials/sl7/prod/external:$PRODUCTS  # Should go away 
setup studio
setup gm2 v9_52_00 -q prof
setup cmake v3_17_1
setup gdb
export USER=lyon
export STUDIO_PROJECT_PATH=/workspaces/irmaData/project
export STUDIO_PROJECT_BUILD_PATH=/workspaces/irmaData/project/build
export STUDIO_PROJECT_SRC_PATH=/workspaces/irmaData/project/src
[[ -z "${FHICL_FILE_PATH}" ]] && FHICL_FILE_PATH=/workspaces/irmaData/project/src || FHICL_FILE_PATH=/workspaces/irmaData/project/src:$FHICL_FILE_PATH
export FHICL_FILE_PATH
LD_LIBRARY_PATH=/workspaces/irmaData/project/build/lib:$LD_LIBRARY_PATH
