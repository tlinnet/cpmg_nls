#!/bin/tcsh
set PLANE=$argv[1]
set METHOD=$argv[2]
set GINI=$argv[3]
set GININI=`echo "$GINI" | cut -f1 -d"_"`
set GINIPLANE=`echo "$GINI" | cut -f2 -d"_"`
set GINIMETHOD=`echo "$GINI" | cut -f3 -d"_"`

set argv[1]=''; set argv[2]=''; set argv[3]=''
set GNINR=( $argv )
set GNINR=( 126_allplanes_coMDD 124_allplanes_coMDD 122_allplanes_coMDD 120_allplanes_coMDD )
set STATS=`awk '{if ( min==""){min=max=$6}; if($6>max) {max=$6}; if($6< min) {min=$6}; total+=$6; count+=1} END {print total/count, min, max, count}' allplanes_${METHOD}.ser`
set STATSAVG=`echo $STATS | cut -f1 -d" "`
set STATSMIN=`echo $STATS | cut -f2 -d" "`
set STATSMAX=`echo $STATS | cut -f3 -d" "`
set STATSNR=`echo $STATS | cut -f4 -d" "`
echo "Avg: $STATSAVG Min: $STATSMIN Max: $STATSMAX Nr: $STATSNR"

cat > corr_coeff.awk <<EOF
#http://people.sc.fsu.edu/~sshanbhag/awk-shell/linreg.awk
#http://www.cyberciti.biz/faq/bash-scripting-using-awk/
#http://en.wikipedia.org/wiki/Simple_linear_regression
#awk -v max=1 '{print \$1, \$2*max}' testdataP225.txt | awk -f corr_coeff.awk
BEGIN{}{ x[NR] = \$1; y[NR] = \$2;  sx += x[NR]; sy += y[NR]; sxx += x[NR]*x[NR]; sxy += x[NR]*y[NR]; syy += y[NR]*y[NR];}
END{ det = NR*sxx - sx*sx; a = (NR*sxy - sx*sy)/det; b = (-sx*sxy+sxx*sy)/det; r = (NR*sxy-sx*sy)/sqrt((NR*sxx-sx*sx)*(NR*syy-sy*sy));
beta=(NR*sxy-sx*sy)/(NR*sxx-sx*sx); alpha=sy/NR-beta/NR*sx; se2=1/(NR*(NR-2))*(NR*syy-sy*sy-beta*beta*(NR*sxx-sx*sx)); sb2=(NR*se2)/(NR*sxx-sx*sx); sa2=sb2/NR*sxx;
for(xi=1;xi<=NR;xi++){xsumsq+=((x[xi]-sx/NR)**2)} ; xvar=xsumsq/(NR); xstdev=sqrt(xvar);
for(yi=1;yi<=NR;yi++){ysumsq+=((y[yi]-sy/NR)**2)} ; yvar=ysumsq/(NR); ystdev=sqrt(xvar);
for(xi=1;xi<=NR;xi++){xysumsq+=((x[xi]-sx/NR)*(y[xi]-sy/NR))} ; xycovar=xysumsq/(NR)
print "NR="NR, "xvar:"xvar, "xstdev="xstdev, "yvar:"yvar, "ystdev="ystdev, "xycovar:"xycovar, "r="r, "r2="r*r, "a=" a, "b="b
#print "sx="sx, "sy="sy, "sxx="sxx, "sxy="sxy, "syy="syy, "beta="beta, "alpha="alpha;
#print "se2="se2, "se="sqrt(se2), "sb2="sb2, "sb="sqrt(sb2), "sa2="sa2, "sa="sqrt(sa2)
}
EOF



set i=9
foreach GNI ( $GNINR )
set GNINI=`echo "$GNI" | cut -f1 -d"_"`
set GNIPLANE=`echo "$GNI" | cut -f2 -d"_"`
set GNIMETHOD=`echo "$GNI" | cut -f3 -d"_"`
echo $GNINI $GNIPLANE $GNIMETHOD 

set METHODSTATS=`awk -v max=$STATSMAX '{print $6*$8/max, $6/max*$'"${i}"'}' ${PLANE}_${METHOD}.ser | awk -f corr_coeff.awk`
echo "$GNINI" "$METHODSTATS" >> ${PLANE}_${METHOD}.stats
@ i++
end
