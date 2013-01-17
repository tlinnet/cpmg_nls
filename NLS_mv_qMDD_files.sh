#!/bin/tcsh

set NCYCPLANESEND=`expr $NCYCPLANES - 1`
#set NCYCPLANESEND=1
set yN=`expr $NINLS \* 2`

foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
mkdir -p ${PLANE}plane.proc
cd ${PLANE}plane.proc
cp -f ../proc.sh .
cp -f ../fidSP.com .
cp -f ../recFT.com .
cp -f ../${PLANE}plane.fid/nls.in .
cp -f ../${PLANE}plane.fid/nls.hdr_3 .
sed -i "/CEXP /{s/yn/nn/}" nls.in
sed -i "s/setenv FID.*/setenv FID ..\/${PLANE}plane/" proc.sh
sed -i "s/setenv NUS_POINTS.*/setenv NUS_POINTS           ${NINLS}/" proc.sh
sed -i "/setenv proc_out /{s/test.dat/test.ft2/}" proc.sh
sed -i "s/-yN.*/-yN          $yN \\/" fidSP.com
sed -i "s/-yT.*/-yT          $NINLS \\/" fidSP.com
#Add CS/MDD parameters
sed  -i  '2{s/^/#Normal parameters\n/}' proc.sh
sed  -i  '2{s/^/setenv CS_niter             10\n/}' proc.sh
sed  -i  '2{s/^/setenv CS_lambda            1.0\n/}' proc.sh
sed  -i  '2{s/^/setenv CS_norm              0\n/}' proc.sh
sed  -i  '2{s/^/setenv CS_alg               IRLS\n/}' proc.sh
sed  -i  '2{s/^/#CS related parameters\n/}' proc.sh
sed  -i  '2{s/^/setenv MDD_NOISE            0.7\n/}' proc.sh
sed  -i  '2{s/^/setenv NITER                50\n/}' proc.sh
sed  -i  '2{s/^/setenv NCOMP                25\n/}' proc.sh
sed  -i  '2{s/^/#MDD related parameters\n/}' proc.sh
#FT
./proc.sh
mv -f test.ft2 ../ft2_data/${NINLS}_${PLANE}_FT.ft2
echo "Done with: ${NINLS}_${PLANE}_FT.ft2"
#CS
sed -i "/setenv METHOD/{s/FT/CS/}" proc.sh
./proc.sh
mv -f test.ft2 ../ft2_data/${NINLS}_${PLANE}_CS.ft2
echo "Done with: ${NINLS}_${PLANE}_CS.ft2"
#MDD
sed -i "/setenv METHOD/{s/CS/MDD/}" proc.sh
./proc.sh
mv -f test.ft2 ../ft2_data/${NINLS}_${PLANE}_MDD.ft2
echo "Done with: ${NINLS}_${PLANE}_MDD.ft2"
#coMDD make ready
sed -i "/set ecode/d" proc.sh
cp proc.sh proc1.sh
sed -i "s/mddnmr4pipeN.sh.*/mddnmr4pipeN.sh 1 23/" proc.sh
sed -i "s/mddnmr4pipeN.sh.*/mddnmr4pipeN.sh 4 5/" proc1.sh
#sed -i "/mddnmr4pipeN.sh /{s/1 2 3 4 5/1 23/}" proc.sh
#sed -i "/mddnmr4pipeN.sh /{s/1 2 3 4 5/4 5/}" proc1.sh
#move nls.in and nls.hdr_3
cp -f nls.in ../ft2_data/${NINLS}_${PLANE}_nls.in
cp -f nls.hdr_3 ../ft2_data/${NINLS}_${PLANE}_nls.hdr_3
cd ..
end

# make coMDD.hd file
rm -rf coMDD
mkdir coMDD
cd coMDD
foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
echo "${PLANE} ${PLANE}plane 1.0 ../${PLANE}plane.proc/MDD/region%02d.mddH ../${PLANE}plane.proc/nls.in 1 2" >> coMDD.hd
end
cd ..
# run coMDD
echo "Running coMDD in $PWD"
set COMDDDIRS=`ls -dtr *.proc`
echo "Parsing coMDD dirs: $COMDDDIRS"
#NLS_coproc.sh c $COMDDDIRS >>&coMDD.log
NLS_coproc.sh c $COMDDDIRS

# Move coMDD results
foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
cd ${PLANE}plane.proc
mv -f test.ft2 ../ft2_data/${NINLS}_${PLANE}_coMDD.ft2
cd ..
end
