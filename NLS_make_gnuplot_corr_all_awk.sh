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

#cat > linreg.awk <<EOF

#EOF

set i=9
foreach GNI ( $GNINR )
set GNINI=`echo "$GNI" | cut -f1 -d"_"`
set GNIPLANE=`echo "$GNI" | cut -f2 -d"_"`
set GNIMETHOD=`echo "$GNI" | cut -f3 -d"_"`
echo $GNINI $GNIPLANE $GNIMETHOD 

#cat > ${PLANE}_${GNINI}_${METHOD}.plt <<EOF

@ i++
end
