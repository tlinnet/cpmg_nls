#!/bin/tcsh

set MODEPAR=`awk '/^array /{f=1;next}f{print $2;exit}' $PROCPARORI | cut -d',' -f1`
if ( $MODEPAR == '"phase' ) then
set MODE=1
endif
if ( $MODEPAR == '"ncyc' ) then
set MODE=0
endif
sort_pseudo3D -in fid -plane $NCYCPLANES -mode $MODE -ni $NI -np $NP
