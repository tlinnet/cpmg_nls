#!/bin/tcsh

# /sbinlab2/tlinnet/Desktop/Link to CPMG_data_test/kte/080716_cpmgDisp_HEWLpH65.fid
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65.fid"
setenv T2 0.1
setenv NIINCR 2

# PROCESS=init for initialization and stopping at 0plane.fid and phase, or PROCESS=run_all
#setenv PROCESS init
setenv PROCESS run_all
NLS_init.sh
