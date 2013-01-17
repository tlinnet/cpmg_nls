#!/bin/tcsh

setenv CPMGFID "/sbinlab2/tlinnet/CPMG_data_test/kte/080716_cpmgDisp_HEWLpH65.fid"
cd $CPMGFID
rm -rf *plane.fid
rm -rf *plane.proc
rm -rf *.fid
rm -rf coMDD
rm proc.sh recFT.com
mv -f procpar_ori procpar
rm -rf nls.in

