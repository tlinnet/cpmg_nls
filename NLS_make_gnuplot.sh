#!/bin/tcsh
set PLANE=$argv[1]
set METHOD=$argv[2]
set GINI=$argv[3]
set GININI=`echo "$GINI" | cut -f1 -d"_"`
set GINIPLANE=`echo "$GINI" | cut -f2 -d"_"`
set GINIMETHOD=`echo "$GINI" | cut -f3 -d"_"`
set argv[1]=''; set argv[2]=''; set argv[3]=''
set GNINR=( $argv )

cat > ${PLANE}_${METHOD}.plt <<EOF
set term postscript eps enhanced color "Helvetica" 14
set title "Relative intensity to ${GININI} ${GINIPLANE} ${GINIMETHOD}\nFor method ${METHOD}, NCYC plane $PLANE"
set size ratio 0.618
set xlabel "Resi number"
set ylabel "Intensity proportionality\nni/ni_{max}"
#set ytics nomirror
#set y2tics border
#set y2label "F-test"
set nokey
#set key autotitle columnheader
set key outside right
#set xrange [0.1:MAXX]
##set xrange[0.5:NRLINES-0.5]
##set xtics 1
set yrange[0:5]
set xtics border in scale 1,0.5 nomirror rotate by -90 font "Helvetica,5"
set output "${PLANE}_${METHOD}.eps"
plot "<(sed -n '14,\${p}' ${PLANE}_${METHOD}.ser)" using 1:(1*\$8):xtic(7) title "$GININI $GINIPLANE $GINIMETHOD",\
EOF
set i=9
foreach GNI ( $GNINR )
    set GNINI=`echo "$GNI" | cut -f1 -d"_"`
    echo '"<(sed -n' "'" '14,${p}' "'" "${PLANE}_${METHOD}.ser" ')"' "using 1:(" '1*$' "${i}):xtic(7) title '$GNINI' ,\\" >> ${PLANE}_${METHOD}.plt
    @ i++
end
set LAST=`rev ${PLANE}_${METHOD}.plt | sed -n '$p' | cut -c3- | rev`
sed -i '$d' ${PLANE}_${METHOD}.plt
echo "$LAST">> ${PLANE}_${METHOD}.plt

gnuplot ${PLANE}_${METHOD}.plt
echo "Done with  ${PLANE}_${METHOD}.eps"
echo "Converting to png: eps2png -resolution 600 ${PLANE}_${METHOD}.eps"
eps2png -resolution 600 ${PLANE}_${METHOD}.eps

