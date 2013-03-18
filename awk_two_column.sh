#!/bin/tcsh -f
#http://people.sc.fsu.edu/~sshanbhag/awk-shell/linreg.awk
#http://www.cyberciti.biz/faq/bash-scripting-using-awk/
#http://en.wikipedia.org/wiki/Simple_linear_regression

set t='"'
if ($#argv < 5) then
    echo "Make statistics from a two column file"
    echo "Usage: $0 file" "'"'$X'"'" "'"'$Y'"' XRATIO YRATIO"
    echo "awk_two_column.sh testdata_wiki.txt" "'"'$1'"'" "'"'$2'"' XRATIO YRATIO 'all'"
    goto done
endif

set FILEIN=$1
#set XCOL=$2
#set YCOL=$3
set XRATIO=$4
set YRATIO=$5
if ( "x$6" == "x" ) then
    set PSTATS='all'
    #echo "PSTATS does not exist. Setting to: $PSTATS"
else
    set PSTATS="$6"
    #echo "Printing for stats: $PSTATS"
endif

if (! -e $FILEIN) then
    echo "$FILEIN does not exist" ; exit 1
endif

#echo "Running: $0 $FILEIN $XCOL $YCOL $XRATIO $YRATIO"

awk "BEGIN{}{x[NR] = ${2}*${XRATIO}; y[NR] = ${3}*${YRATIO}; res[NR]=y[NR]-x[NR] \
if(minx==${t}${t}){minx=maxx=x[NR];miny=maxy=y[NR];minres=maxres=res[NR];minxi=maxxi=minyi=maxyi=minresi=maxresi=1}; \
if(x[NR]<minx){minx=x[NR];minxi=NR};if(x[NR]>maxx){maxx=x[NR];maxxi=NR}; \
if(y[NR]<miny){miny=y[NR];minyi=NR};if(y[NR]>maxy){maxy=y[NR];maxyi=NR}; \
if(res[NR]<minres){minres=res[NR];minresi=NR};if(res[NR]>maxres){maxres=res[NR];maxresi=NR}; \
sx += x[NR]; sy += y[NR]; sxx += x[NR]*x[NR]; syy += y[NR]*y[NR]; sres += res[NR]; sxy += x[NR]*y[NR]} \
END{\
det = NR*sxx - sx*sx; a = (NR*sxy - sx*sy)/det; b = (-sx*sxy+sxx*sy)/det; r = (NR*sxy-sx*sy)/sqrt((NR*sxx-sx*sx)*(NR*syy-sy*sy)); \
beta=(NR*sxy-sx*sy)/(NR*sxx-sx*sx); alpha=sy/NR-beta/NR*sx; se2=1/(NR*(NR-2))*(NR*syy-sy*sy-beta*beta*(NR*sxx-sx*sx)); sb2=(NR*se2)/(NR*sxx-sx*sx); sa2=sb2/NR*sxx; \
for(xi=1;xi<=NR;xi++){\
xsumsq+=((x[xi]-sx/NR)**2);\
ysumsq+=((y[xi]-sy/NR)**2);\
ressumsq+=((y[xi]-x[xi]-sres/NR)**2)\
xysumsq+=((x[xi]-sx/NR)*(y[xi]-sy/NR))\
} \
xvar=xsumsq/(NR); xstdev=sqrt(xvar); \
yvar=ysumsq/(NR); ystdev=sqrt(yvar); \
resvar=ressumsq/(NR); resstdev=sqrt(resvar); \
xycovar=xysumsq/(NR); \
print ${t}NR= $t NR, ${t}minx= $t minx, ${t}minxi= $t minxi, ${t}maxx= $t maxx, ${t}maxxi= $t maxxi, ${t}miny= $t miny, ${t}minyi= $t minyi, ${t}maxy= $t maxy, ${t}maxyi= $t maxyi, ${t}minres= $t minres, ${t}minresi= $t minresi, ${t}maxres= $t maxres, ${t}maxresi= $t maxresi \
print ${t}xvar= $t xvar, ${t}xstdev= $t xstdev, ${t}yvar= $t yvar, ${t}ystdev= $t ystdev, ${t}resvar= $t resvar, ${t}resstdev= $t resstdev, ${t}xycovar= $t xycovar \
print ${t}r= $t r, ${t}r2= $t r*r, ${t}a= $t a, ${t}b= $t b, ${t}alpha= $t alpha, ${t}beta= $t beta \
print ${t}sx= $t sx, ${t}sy= $t sy, ${t}sxx= $t sxx, ${t}sxy= $t sxy, ${t}syy= $t syy \
print ${t}se2= $t se2, ${t}sb2= $t sb2, ${t}sa2= $t sa2, ${t}sresNR= $t sres/NR \
}" $FILEIN > ${FILEIN}.stats
set STATS=`cat ${FILEIN}.stats`
#rm ${FILEIN}.stats

set NRSTATS=${#STATS}
set NRSTATSEND=`expr $NRSTATS - 1`

foreach i (`seq 1 2 $NRSTATSEND`)
    set j=`expr $i + 1`
    if ( "$PSTATS" == "all") then
        echo -n "$STATS[$i]$i $STATS[$j] "
    else
        foreach PSTAT ( $PSTATS )
            if ( $STATS[$i] == ${PSTAT}=) then
                echo -n "$STATS[$i] $STATS[$j] "
            else if ( $STATS[$i] == ${PSTAT}) then
                echo -n "$STATS[$i] $STATS[$j] "
            endif
        end
    endif    
end
echo

done: ; exit 0
error: ; exit 1
