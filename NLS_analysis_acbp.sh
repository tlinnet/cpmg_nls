#!/bin/tcsh

setenv SPARKYPEAKLIST peaks.list
############### cpmgDisp_HEWLpH65 #############################
setenv CPMGFID "/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_normal.fid"
#NLS_doanalysis.sh

#Running for CS30 and MDD500
setenv CPMGFID "/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS30_MDD500.fid"
#NLS_doanalysis.sh

#Running for CS100 and MDD500
setenv CPMGFID "/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS100_MDD500.fid"
NLS_doanalysis.sh

###############
