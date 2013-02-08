#!/bin/tcsh

set EXPRANGE=( $argv )
set yN=`expr $NINLS \* 2`

foreach EXP ( $EXPRANGE )
set EXPDIR=`echo $EXP | cut -d'.' -f1`
mkdir -p ${EXPDIR}.proc
cd ${EXPDIR}.proc
cp -f ../proc.sh .
cp -f ../fidSP.com .
cp -f ../recFT.com .
cp -f ../${EXPDIR}.fid/nls.in .
cp -f ../${EXPDIR}.fid/nls.hdr_3 .
sed -i "/CEXP /{s/yn/nn/}" nls.in
sed -i "s/setenv FID.*/setenv FID ..\/${EXPDIR}/" proc.sh
sed -i "s/setenv NUS_POINTS.*/setenv NUS_POINTS           ${NINLS}/" proc.sh
sed -i "/setenv proc_out /{s/test.dat/test.ft2/}" proc.sh
sed -i "s/-yN.*/-yN          $yN \\/" fidSP.com
sed -i "s/-yT.*/-yT          $NINLS \\/" fidSP.com
#Add CS/MDD parameters
sed  -i  "2{s/^/#Normal parameters\n/}" proc.sh
sed  -i  "2{s/^/setenv CS_niter             ${CS_niter}\n/}" proc.sh
sed  -i  "2{s/^/setenv CS_lambda            ${CS_lambda}\n/}" proc.sh
sed  -i  "2{s/^/setenv CS_norm              ${CS_norm}\n/}" proc.sh
sed  -i  "2{s/^/setenv CS_alg               ${CS_alg}\n/}" proc.sh
sed  -i  "2{s/^/#CS related parameters\n/}" proc.sh
sed  -i  "2{s/^/setenv MDD_NOISE            ${MDD_NOISE}\n/}" proc.sh
sed  -i  "2{s/^/setenv NITER                ${NITER}\n/}" proc.sh
sed  -i  "2{s/^/setenv NCOMP                ${NCOMP}\n/}" proc.sh
sed  -i  "2{s/^/#MDD related parameters\n/}" proc.sh
#FT
./proc.sh
mv -f test.ft2 ../$FTDATA/${NINLS}_${EXPDIR}_FT.ft2
echo "Done with: ${NINLS}_${EXPDIR}_FT.ft2"
#CS
sed -i "/setenv METHOD/{s/FT/CS/}" proc.sh
./proc.sh
mv -f test.ft2 ../$FTDATA/${NINLS}_${EXPDIR}_CS.ft2
echo "Done with: ${NINLS}_${EXPDIR}_CS.ft2"
#MDD
sed -i "/setenv METHOD/{s/CS/MDD/}" proc.sh
./proc.sh
mv -f test.ft2 ../$FTDATA/${NINLS}_${EXPDIR}_MDD.ft2
echo "Done with: ${NINLS}_${EXPDIR}_MDD.ft2"
#coMDD make ready
sed -i "/set ecode/d" proc.sh
cp proc.sh proc1.sh
sed -i "s/mddnmr4pipeN.sh.*/mddnmr4pipeN.sh 1 23/" proc.sh
sed -i "s/mddnmr4pipeN.sh.*/mddnmr4pipeN.sh 4 5/" proc1.sh
#move nls.in and nls.hdr_3
cp -f nls.in ../$FTDATA/${NINLS}_${EXPDIR}_nls.in
cp -f nls.hdr_3 ../$FTDATA/${NINLS}_${EXPDIR}_nls.hdr_3
cd ..
end

# make coMDD.hd file
rm -rf coMDD
mkdir coMDD
cd coMDD
foreach EXP ( $EXPRANGE )
set EXPDIR=`echo $EXP | cut -d'.' -f1`
echo "${EXPDIR} ${EXPDIR} 1.0 ../${EXPDIR}.proc/MDD/region%02d.mddH ../${EXPDIR}.proc/nls.in 1 2" >> coMDD.hd
end
cd ..
# run coMDD
echo "Running coMDD in $PWD"
set COMDDDIRS=`ls -dtr *.proc`
echo "Parsing coMDD dirs: $COMDDDIRS"
NLS_coproc.sh c $COMDDDIRS

# Move coMDD results
foreach EXP ( $EXPRANGE )
set EXPDIR=`echo $EXP | cut -d'.' -f1`
cd ${EXPDIR}.proc
mv -f test.ft2 ../$FTDATA/${NINLS}_${EXPDIR}_coMDD.ft2
cd ..
end
