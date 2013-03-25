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

# PROCESS=init for initialization and stopping at 0plane.fid and phase, or PROCESS=run_all
############## SOD1 cpmg_disp_sod1d90a_060518 #############################
setenv CPMGFID "/sbinlab2/tlinnet/Desktop/haddock_tlinnet/kte/SOD1/cpmg_disp_sod1d90a_060518/cpmg_disp_sod1d90a_060518_normal.fid"
#NLS_init.sh init
NLS_init.sh run_all >> $CPMGFID/logfile.log

############## SOD1 cpmg_disp_sod1d90a_060521 #############################
setenv CPMGFID "/sbinlab2/tlinnet/Desktop/haddock_tlinnet/kte/SOD1/cpmg_disp_sod1d90a_060521/cpmg_disp_sod1d90a_060521_normal.fid"
#NLS_init.sh init
NLS_init.sh run_all >> $CPMGFID/logfile.log

############## SOD1 cpmg_disp_sod1WT_050630 #############################
setenv CPMGFID "/sbinlab2/tlinnet/Desktop/haddock_tlinnet/kte/SOD1/cpmg_disp_sod1WT_050630/cpmg_disp_sod1WT_050630_normal.fid"
#NLS_init.sh init
NLS_init.sh run_all >> $CPMGFID/logfile.log

############## SOD1 cpmg_disp_sod1WT_050720 #############################
setenv CPMGFID "/sbinlab2/tlinnet/Desktop/haddock_tlinnet/kte/SOD1/cpmg_disp_sod1WT_050720/cpmg_disp_sod1WT_050720_normal.fid"
#NLS_init.sh init
NLS_init.sh run_all >> $CPMGFID/logfile.log
