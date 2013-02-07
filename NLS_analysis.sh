#!/bin/tcsh

############### cpmgDisp_HEWLpH65 #############################
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal.fid"
NLS_doanalysis.sh
#Running for NoRandom schedule
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_NoRandom.fid"
NLS_doanalysis.sh
#Running for CS30 and MDD500
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500.fid"
NLS_doanalysis.sh

#Running for CS30 and MDD500 for NoRandom schedule
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500_NoRandom.fid"
NLS_doanalysis.sh

###############
