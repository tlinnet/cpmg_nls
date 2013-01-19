#!/bin/tcsh
set PLANE=$argv[1]
set METHOD=$argv[2]
set GINI=$argv[3]
set GININI=`echo "$GINI" | cut -f1 -d"_"`
set GINIPLANE=`echo "$GINI" | cut -f2 -d"_"`
set GINIMETHOD=`echo "$GINI" | cut -f3 -d"_"`

set argv[1]=''; set argv[2]=''; set argv[3]=''
set GNINR=( $argv )
#set GNINR=( 126_0_coMDD 124_0_coMDD 122_0_coMDD 120_0_coMDD 120_0_coMDD )
set STATS=`sed -n '14,${p}' ${PLANE}_${METHOD}.ser | awk '{if ( min==""){min=max=$6}; if($6>max) {max=$6}; if($6< min) {min=$6}; total+=$6; count+=1} END {print total/count, min, max, count}'`
set STATSAVG=`echo $STATS | cut -f1 -d" "`
set STATSMIN=`echo $STATS | cut -f2 -d" "`
set STATSMAX=`echo $STATS | cut -f3 -d" "`
set STATSNR=`echo $STATS | cut -f4 -d" "`

set i=9
foreach GNI ( $GNINR )
set GNINI=`echo "$GNI" | cut -f1 -d"_"`
set GNIPLANE=`echo "$GNI" | cut -f2 -d"_"`
set GNIMETHOD=`echo "$GNI" | cut -f3 -d"_"`


cat > ${PLANE}_${GNINI}_${METHOD}.plt <<EOF
set term postscript eps enhanced color "Helvetica" 14
set title "Intensity correlation plot for method ${METHOD}\nIntensity normalized to maximum of ${GININI} ${GINIPLANE} ${GINIMETHOD}"
set size ratio 0.618
set xlabel "Intensity ni=${GININI} NCYC plane=${GINIPLANE} for ${GINIMETHOD}"
set ylabel "Intensity ni=${GNINI} NCYC plane=${GNIPLANE} for ${GNIMETHOD}"
#set ytics nomirror
#set y2tics border
#set y2label "F-test"
set nokey
#set key autotitle columnheader
#set key outside right
set xrange [0:]
set yrange[0:]
#set xtics border in scale 1,0.5 mirror rotate by -90 font "Helvetica,5"
#set ytics border in scale 1,0.5 mirror font "Helvetica,5"
set output "${PLANE}_${GNINI}_${METHOD}.eps"

set fit errorvariables
f(x)=a*x + b
fit f(x) "<(sed -n '14,\${p}' ${PLANE}_${METHOD}.ser)" using (\$6*\$8 /${STATSMAX}):(\$6*\$$i /${STATSMAX}) via a, b
f_FIT_NDF=FIT_NDF
f_FIT_STDFIT=FIT_STDFIT
f_FIT_WSSR=FIT_WSSR
f_WSSR_NDF=FIT_WSSR/FIT_NDF

set label sprintf("f(x)=a*x") at graph 0.8,0.80 font "Helvetica,9"
set label sprintf("a=%1.4f +- %1.5f",a,a_err) at graph 0.8,0.75 font "Helvetica,9"
set label sprintf("b=%1.4f +- %1.5f",b,b_err) at graph 0.8,0.70 font "Helvetica,9"
set label sprintf("NDF=%1.4f",f_FIT_NDF) at graph 0.8,0.65 font "Helvetica,9"
set label sprintf("STDFIT=%1.7f",f_FIT_STDFIT) at graph 0.8,0.60 font "Helvetica,9"
set label sprintf("WSSR=%1.7f",f_FIT_WSSR) at graph 0.8,0.55 font "Helvetica,9"
set label sprintf("WSSR/NDF=%1.7f",f_WSSR_NDF) at graph 0.8,0.50 font "Helvetica,9"

plot "<(sed -n '14,\${p}' ${PLANE}_${METHOD}.ser)" using (\$6*\$8 /${STATSMAX}):(\$6*\$$i /${STATSMAX}) title "${GNINI} ${GNIPLANE} ${GNIMETHOD}",\
f(x) title "Correlation"
EOF

gnuplot ${PLANE}_${GNINI}_${METHOD}.plt
echo "Done with  ${PLANE}_${GNINI}_${METHOD}.eps"
echo "Converting to png: eps2png -resolution 400 ${PLANE}_${GNINI}_${METHOD}.eps"
eps2png -resolution 400 ${PLANE}_${GNINI}_${METHOD}.eps

@ i++
end

