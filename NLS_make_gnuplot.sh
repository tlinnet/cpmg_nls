#!/bin/tcsh
set PLANE=$argv[1]
set METHOD=$argv[2]
set GINI=$argv[3]
set argv[1]=''; set argv[2]=''; set argv[3]=''
set GNINR=( $argv )

cat > ${PLANE}_${METHOD}.plt <<EOF
set term postscript eps enhanced color "Helvetica" 14
set title "For method ${METHOD}"
set size ratio 0.3
set xlabel "Resi number"
set ylabel "Relative intensity to ${NI}_${PLANE}_FT"
#set ytics nomirror
#set y2tics border
#set y2label "F-test"
set nokey
#set key autotitle columnheader
#set key outside right
#set xrange [0.1:MAXX]
##set xrange[0.5:NRLINES-0.5]
##set xtics 1
#set yrange[0:0.025]
set xtics border in scale 1,0.5 nomirror rotate by -90 font "Helvetica,5"
set output "${PLANE}_${METHOD}.eps"
plot "<(sed -n '14,\${p}' ${PLANE}_${METHOD}.ser)" using 1:(1*\$8):xtic(7) title "$GINI",\
EOF
set i=8
foreach GNI ( $GNINR )
    @ i++
    echo '"<(sed -n' "'" '14,${p}' "'" "${PLANE}_${METHOD}.ser" ')"' "using 1:(" '1*$' "${i}):xtic(7) title '$GNI' ,\\" >> ${PLANE}_${METHOD}.plt
end
set LAST=`rev ${PLANE}_${METHOD}.plt | sed -n '$p' | cut -c3- | rev`
sed -i '$d' ${PLANE}_${METHOD}.plt
echo "$LAST">> ${PLANE}_${METHOD}.plt

gnuplot ${PLANE}_${METHOD}.plt
echo "Converting to png: eps2png -resolution 200 ${PLANE}_${METHOD}.eps"
eps2png -resolution 200 ${PLANE}_${METHOD}.eps

