#!/bin/tcsh

if ( $HOST != "haddock") then
echo "You have to run on haddock"
exit
endif

cd $TRHOFIDS
set EXPFILES=`ls -vd $TRHOFIDS/*.fid`
setenv EXPFILESNR ${#EXPFILES}
#setenv EXPFILESNR 2
set EXPRANGE = ""
foreach I (`seq 1 1 $EXPFILESNR`)
    set EXPTEMP=`basename "$EXPFILES[$I]"`
    set EXPRANGE="${EXPRANGE} $EXPTEMP "
end

setenv PROCESS $1
setenv PROCPAR procpar
setenv PROCPARORI ${PROCPAR}_ori
setenv EXPINI $EXPFILES[1]
echo "expini: $EXPINI"
cp -pn $EXPINI/$PROCPAR $TRHOFIDS/$PROCPAR
setenv NI `awk '/^ni /{f=1;next}f{print $2;exit}' $PROCPAR`
setenv NIEND `expr $NI - 1`
setenv NP `awk '/^np /{f=1;next}f{print $2;exit}' $PROCPAR`
setenv sw1 `sed -n "/sw1 /{n;p}" $PROCPAR | awk '{print $2}'`
setenv sw `sed -n "/sw /{n;p}" $PROCPAR | awk '{print $2}'`
setenv FTDATA ft2_data

if ( $PROCESS == "init") then
cp -pn $PROCPAR $PROCPARORI
setenv NINLS $NI
NLS_init_nls.sh
TRHO_mv_truncate_fid.sh $EXPRANGE

qMDD22 $EXPINI  
cat << EOF
##########################
1) Now phase your spectrum for the 0plane.
In fidSP.com:
| nmrPipe -fn PS -p0 0 -p1 0 -di  \

2) Correct the Fourier Transform (FT) if the spectrum is mirrored in the horizontal plane. 
In recFT.com:
| nmrPipe  -fn FT -neg \


Then rerun your NLS_proc.sh with PROCESS=run_all

##########################
EOF

else if ( $PROCESS == "run_all") then
set EXPINIDIR=`echo $EXPINI | cut -d'.' -f1`
mv -n ${EXPINIDIR}.proc/proc.sh . ; mv -n ${EXPINIDIR}.proc/fidSP.com . ; mv -n ${EXPINIDIR}.proc/recFT.com .
mkdir -p $FTDATA

foreach NLSNI (`seq 0 $NIINCR $NIEND`)
setenv NINLS `expr $NI - $NLSNI`
echo "Processing NI=${NINLS}"
TRHO_mv_truncate_fid.sh $EXPRANGE
TRHO_mv_qMDD_files.sh $EXPRANGE
end


endif
