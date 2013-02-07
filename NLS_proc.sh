#!/bin/tcsh

# Parameters passed to processing
#CS related parameters
setenv CS_alg IRLS
setenv CS_norm 0
setenv CS_lambda 1.0
setenv CS_niter 10
#MDD related parameters
setenv NCOMP 25
setenv NITER 50
setenv MDD_NOISE 0.7

# Parameters passed to nussampler
setenv T2 0.1
setenv NIINCR 2

# Parameters for data directory
# /sbinlab2/tlinnet/Desktop/Link to CPMG_data_test/kte/080716_cpmgDisp_HEWLpH65.fid
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65.fid"
setenv NUSRANDOM "y"

# PROCESS=init for initialization and stopping at 0plane.fid and phase, or PROCESS=run_all
NLS_init.sh init
#setenv PROCESS run_all >> $CPMGFID/logfile.log

