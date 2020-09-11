kx509
voms-proxy-init -rfc -noregen -voms=fermilab:/fermilab/gm2/Role=Analysis
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
setup fife_utils
export SAM_EXPERIMENT=gm2
export EXPERIMENT=gm2