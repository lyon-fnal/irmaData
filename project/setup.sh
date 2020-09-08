#!/usr/bin/env bash
source /cvmfs/gm2.opensciencegrid.org/prod/g-2/setup
export PRODUCTS=${WORKSPACE_FOLDER}/project/products:/cvmfs/gm2.opensciencegrid.org/specials/sl7/prod/g-2:/cvmfs/gm2.opensciencegrid.org/specials/sl7/prod/external:$PRODUCTS  # Should go away 
setup studio
setup gm2 v9_52_00 -q prof
setup cmake v3_17_1
setup hep_hpc v0_13_01 -q e15
setup gdb
export STUDIO_PROJECT_PATH=${WORKSPACE_FOLDER}/project
export STUDIO_PROJECT_BUILD_PATH=${WORKSPACE_FOLDER}/project/build
export STUDIO_PROJECT_SRC_PATH=${WORKSPACE_FOLDER}/project/src
[[ -z "${FHICL_FILE_PATH}" ]] && FHICL_FILE_PATH=${WORKSPACE_FOLDER}/project/src || FHICL_FILE_PATH=${WORKSPACE_FOLDER}/project/src:$FHICL_FILE_PATH
export FHICL_FILE_PATH
LD_LIBRARY_PATH=${WORKSPACE_FOLDER}/project/build/lib:$LD_LIBRARY_PATH
