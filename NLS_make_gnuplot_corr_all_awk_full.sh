#!/bin/tcsh
set PLANE=$argv[1]
set METHOD=$argv[2]
set GINI=$argv[3]
set GININI=`echo "$GINI" | cut -f1 -d"_"`
set GINIPLANE=`echo "$GINI" | cut -f2 -d"_"`
set GINIMETHOD=`echo "$GINI" | cut -f3 -d"_"`
set t="'"

set argv[1]=''; set argv[2]=''; set argv[3]=''
set GNINR=( $argv )
#set GNINR=( 126_allplanes_coMDD 124_allplanes_coMDD 122_allplanes_coMDD 120_allplanes_coMDD )
set STATS=`awk_two_column.sh ${PLANE}_${METHOD}.ser '$6*$8' '$6*$9' 1 1 'maxx a'`
set STATSMAX=`echo $STATS | cut -f2 -d" "`
set STATSA=`echo $STATS | cut -f4 -d" "`

set i=9
foreach GNI ( $GNINR )
set GNINI=`echo "$GNI" | cut -f1 -d"_"`
set GNIPLANE=`echo "$GNI" | cut -f2 -d"_"`
set GNIMETHOD=`echo "$GNI" | cut -f3 -d"_"`
echo "Finding stats for:" $GNINI $GNIPLANE $GNIMETHOD 

#echo awk_two_column.sh ${PLANE}_${METHOD}.ser '$6*$8' '$6*$'"${i}" 1/$STATSMAX 1/$STATSMAX/$STATSA 'all'
set METHODSTATS=`awk_two_column.sh ${PLANE}_${METHOD}.ser '$6*$8' '$6*$'"${i}" 1/$STATSMAX 1/$STATSMAX/$STATSA 'all'`
echo "$GNINI" "$METHODSTATS" >> ${PLANE}_${METHOD}.stats
@ i++
end

cat > ${PLANE}_${METHOD}.plt <<EOF
set term postscript eps enhanced color "Helvetica" 14
set title "Statistics plot for method ${METHOD}\nIntensity normalized to maximum of ${GININI} ${GINIPLANE} ${GINIMETHOD}, and a=${STATSA}"
set size ratio 0.618
set xlabel "Number of random ni out of ni_{max}"
set ylabel "Correlation Coefficient R^2\nIntensity proportionality a, f(x)=a*x+b"
set ytics nomirror
set y2tics border
set y2label "Residual stdev {/Symbol s}"
#set y2label "Covariance {/Symbol s}_{xy}"
set nokey
#set key autotitle columnheader
set key bottom center outside font "Helvetica,7"
set xrange [-${GININI}:]
set yrange[-1:1.1]
set y2range[0:0.1]
#set xtics border in scale 1,0.5 mirror rotate by -90 font "Helvetica,5"
#set ytics border in scale 1,0.5 mirror font "Helvetica,5"
set output "${PLANE}_${METHOD}.eps"

#set fit errorvariables
#f(x)=a*x + b
#fit f(x) "${PLANE}_${METHOD}.ser" using (\$6*\$8 /${STATSMAX}):(\$6*\$$i /${STATSMAX}) via a, b
#f_FIT_NDF=FIT_NDF
#f_FIT_STDFIT=FIT_STDFIT
#f_FIT_WSSR=FIT_WSSR
#f_WSSR_NDF=FIT_WSSR/FIT_NDF

#set label sprintf("f(x)=a*x") at graph 0.8,0.80 font "Helvetica,9"
#set label sprintf("a=%1.2f +- %1.4f",a,a_err) at graph 0.8,0.70 font "Helvetica,9"
#set label sprintf("b=%1.2f +- %1.4f",b,b_err) at graph 0.8,0.60 font "Helvetica,9"
#set label sprintf("NDF=%1.2f",f_FIT_NDF) at graph 0.8,0.50 font "Helvetica,9"
#set label sprintf("STDFIT=%1.2f",f_FIT_STDFIT) at graph 0.8,0.40 font "Helvetica,9"
#set label sprintf("WSSR=%1.2f",f_FIT_WSSR) at graph 0.8,0.30 font "Helvetica,9"
#set label sprintf("WSSR/NDF=%1.2f",f_WSSR_NDF) at graph 0.8,0.20 font "Helvetica,9"

plot "${PLANE}_${METHOD}.stats" using (-1*\$1):(\$71):(\$21):(\$25) title "Residual average" with yerrorbars,\
"${PLANE}_${METHOD}.stats" using (-1*\$1):(\$47) title "Intensity proportionality a, f(x)=a*x+b",\
"${PLANE}_${METHOD}.stats" using (-1*\$1):(\$45) title "Correlation Coefficient R^2",\
"${PLANE}_${METHOD}.stats" using (-1*\$1):(\$39) title "Y2: Residual stdev {/Symbol s}" axis x1y2,\
"${PLANE}_${METHOD}.stats" using (-1*\$1):(\$41) title "Y2: Covariance {/Symbol s}_{xy}" axis x1y2
EOF

gnuplot ${PLANE}_${METHOD}.plt
echo "Done with  ${PLANE}_${METHOD}.eps"
echo "Converting to png: eps2png -resolution 400 ${PLANE}_${METHOD}.eps"
eps2png -resolution 400 ${PLANE}_${METHOD}.eps
