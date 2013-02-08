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

############### cpmgDisp_HEWLpH65 #############################
# PROCESS=init for initialization and stopping at 0plane.fid and phase, or PROCESS=run_all
setenv TRHOFIDS "/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_normal"
setenv NUSRANDOM "y"
#TRHO_init.sh init
TRHO_init.sh run_all >> $TRHOFIDS/logfile.log

#Running for NoRandom schedule
setenv TRHOFIDS "/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_NoRandom"
setenv NUSRANDOM "n"
TRHO_init.sh run_all >> $TRHOFIDS/logfile.log
