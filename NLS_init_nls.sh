#!/bin/tcsh

set SPARSE=`sed -n "/SPARSE /{n;p}" $PROCPAR | awk '{print $2}'`
if ( $SPARSE == "" ) then
echo "SPARSE " >> $PROCPAR
echo '1 "y"' >> $PROCPAR
echo '0' >> $PROCPAR
endif
sed -i "/SPARSE /{n; s/n/y/}" $PROCPAR

set f1coef=`sed -n "/f1coef /{n;p}" $PROCPAR | awk '{print $1}'`
if ( $f1coef == "" ) then
echo "f1coef " >> $PROCPAR
echo '1 "1 0 -1 0 0 -1 0 -1"' >> $PROCPAR
echo '0' >> $PROCPAR
endif

cat > nls.in << EOF
NDIM  2
SPARSE y
sptype shuffle
seed   4321
CEXP   yn
CT_SP  nn
NIMAX  $NI  1
NI     $NI  1
SW     $sw1 $sw
T2      $T2 0.05  1
Jsp     0 0 0
phase   0 0 0 
f180    nn
EOF
