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
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal.fid"
setenv NUSRANDOM "n"
#NLS_init.sh init
NLS_init.sh run_all >> $CPMGFID/logfile.log

#Running for NoRandom schedule
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_NoRandom.fid"
setenv NUSRANDOM "n"
NLS_init.sh run_all >> $CPMGFID/logfile.log

#Running for CS30 and MDD500
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500.fid"
setenv NUSRANDOM "y"
setenv CS_niter 30
setenv NITER 500
NLS_init.sh run_all >> $CPMGFID/logfile.log

#Running for CS30 and MDD500 for NoRandom schedule
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500_NoRandom.fid"
setenv NUSRANDOM "n"
setenv CS_niter 30
setenv NITER 500
NLS_init.sh run_all >> $CPMGFID/logfile.log

###############
