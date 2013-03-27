#!/bin/tcsh

setenv SPARKYPEAKLIST peaks.list
############## SOD1 cpmg_disp_sod1d90a_060518 #############################
setenv CPMGFID "/home/tlinnet/kte/SOD1/cpmg_disp_sod1d90a_060518/cpmg_disp_sod1d90a_060518_normal.fid"
NLS_doanalysis.sh

############## SOD1 cpmg_disp_sod1d90a_060521 #############################
setenv CPMGFID "/home/tlinnet/kte/SOD1/cpmg_disp_sod1d90a_060521/cpmg_disp_sod1d90a_060521_normal.fid"
NLS_doanalysis.sh

############## SOD1 cpmg_disp_sod1WT_050630 #############################
setenv CPMGFID "/home/tlinnet/kte/SOD1/cpmg_disp_sod1WT_050630/cpmg_disp_sod1WT_050630_normal.fid"
NLS_doanalysis.sh

############## SOD1 cpmg_disp_sod1WT_050720 #############################
setenv CPMGFID "/home/tlinnet/kte/SOD1/cpmg_disp_sod1WT_050720/cpmg_disp_sod1WT_050720_normal.fid"
NLS_doanalysis.sh
