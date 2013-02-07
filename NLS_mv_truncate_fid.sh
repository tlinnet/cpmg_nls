#!/bin/tcsh

set NCYCPLANESEND=`expr $NCYCPLANES - 1`
#set NCYCPLANESEND=1

foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
mkdir -p ${PLANE}plane.fid
cd ${PLANE}plane.fid
cp -n ../${PLANE}.fid ./fid_ori
cp -f ../nls.in .
cp -f ../procpar .
sed -i "/ni /{n; s/$NI/$NINLS/}" $PROCPAR
sed -i "/array /{n; s/phase,ncyc/phase/}" $PROCPAR
sed -i "/NI /{s/$NI/$NINLS/}" nls.in
if ( $NUSRANDOM == "y") then
set random1=`head -c 1 /dev/urandom | od -t u1 | cut -c9-`
set random2=`head -c 1 /dev/urandom | od -t u1 | cut -c9-`
set random3=`expr $random1 \* $random2` 
sed -i "/seed /{s/4321/$random3/}" nls.in
endif
nussampler nls.in
sort -n nls.hdr_3 > tmp; rm nls.hdr_3; mv tmp nls.hdr_3
makeNUSdata.pl fid_ori nls.hdr_3 fid $NP $NI
cd ..
end
