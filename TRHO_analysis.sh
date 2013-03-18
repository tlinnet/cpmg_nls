#!/bin/tcsh

#setenv SPARKYPEAKLIST sparky.peaks
setenv SPARKYPEAKLIST sparky.list
############### cpmgDisp_HEWLpH65 #############################
#setenv TRHOFIDS "/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_normal"
#TRHO_doanalysis.sh
#Running for 

#Running for CS30 and MDD500
setenv TRHOFIDS "/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS30_MDD500"
TRHO_doanalysis.sh

#Running for CS100 and MDD500
setenv TRHOFIDS "/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS100_MDD500"
TRHO_doanalysis.sh


###############
