#!/bin/tcsh

if ( $HOST != "haddock") then
echo "You have to run on haddock"
exit
endif

cd $CPMGFID
setenv PROCPAR procpar
setenv PROCPARORI ${PROCPAR}_ori
setenv NCYCPLANES `awk '/^ncyc /{f=1;next}f{print $1;exit}' $PROCPAR`
setenv NI `awk '/^ni /{f=1;next}f{print $2;exit}' $PROCPAR`
setenv NIEND `expr $NI - 1`
setenv NP `awk '/^np /{f=1;next}f{print $2;exit}' $PROCPAR`
setenv sw1 `sed -n "/sw1 /{n;p}" $PROCPAR | awk '{print $2}'`
setenv sw `sed -n "/sw /{n;p}" $PROCPAR | awk '{print $2}'`

if ( $PROCESS == "init") then
cp -pn $PROCPAR $PROCPARORI
NLS_init_fid.sh

setenv NINLS $NI
NLS_init_nls.sh
NLS_mv_truncate_fid.sh

#qMDD nogui 0plane.fid --METHOD FT --mddargs FST_PNT_PPM 11 ROISW 5 proc_out test.ft2 MDDTHREADS 22
qMDD 0plane.fid  
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
mv -n 0plane.proc/proc.sh . ; mv -n 0plane.proc/fidSP.com . ; mv -n 0plane.proc/recFT.com .
mkdir -p ft2_data

foreach NLSNI (`seq 0 $NIINCR $NIEND`)
setenv NINLS `expr $NI - $NLSNI`
echo "Processing NI=${NINLS}"
NLS_mv_truncate_fid.sh
NLS_mv_qMDD_files.sh
end

endif
