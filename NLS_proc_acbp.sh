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
setenv NUSRANDOM "y"

############### cpmgDisp_HEWLpH65 #############################
# PROCESS=init for initialization and stopping at 0plane.fid and phase, or PROCESS=run_all
setenv CPMGFID "/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_normal.fid"
#NLS_init.sh init
#NLS_init.sh run_all >> $CPMGFID/logfile.log

#Running for CS30 and MDD500
setenv CPMGFID "/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS30_MDD500.fid"
setenv CS_niter 30
setenv NITER 500
#NLS_init.sh run_all >> $CPMGFID/logfile.log

#Running for CS100 and MDD500
setenv CPMGFID "/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS100_MDD500.fid"
setenv CS_niter 100
setenv NITER 500
#NLS_init.sh run_all >> $CPMGFID/logfile.log

# Running for normal, no Noise
setenv MDD_NOISE 0.0
setenv CS_niter 10
setenv NITER 50
setenv CPMGFID "/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_normal_noNoise.fid"
NLS_init.sh run_all >> $CPMGFID/logfile.log

###############
