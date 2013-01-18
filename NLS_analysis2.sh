#!/bin/tcsh

# /sbinlab2/tlinnet/Desktop/Link to CPMG_data_test/kte/080716_cpmgDisp_HEWLpH65.fid
setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65.fid"
setenv NIINCR 2
setenv ANADIR $CPMGFID/ft2_data/analysis2
setenv SPARKYPEAKLIST sparky.list
cd $CPMGFID

setenv PROCPAR procpar
setenv PROCPARORI ${PROCPAR}_ori
setenv NCYCPLANES `awk '/^ncyc /{f=1;next}f{print $1;exit}' $PROCPARORI`
set NCYCPLANESEND=`expr $NCYCPLANES - 1`
setenv NI `awk '/^ni /{f=1;next}f{print $2;exit}' $PROCPARORI`
setenv NIEND `expr $NI - 1`
setenv NP `awk '/^np /{f=1;next}f{print $2;exit}' $PROCPARORI`

mkdir -p $ANADIR
cd $ANADIR
echo "working in $PWD"
cp -n $CPMGFID/$SPARKYPEAKLIST .
echo "Making peaks.dat"
echo "NLS_stPeakList.pl ../128_0_FT.ft2 $SPARKYPEAKLIST > peaks.dat\n"
NLS_stPeakList.pl ../128_0_FT.ft2 $SPARKYPEAKLIST > peaks.dat


#set PROCESS="int_resi"
#set PROCESS="int_corr_ft_method"
set PROCESS="int_corr_ft_method_all"

if ( $PROCESS == "int_resi") then
#set NCYCPLANESEND=0
set METHODS = "FT CS MDD coMDD"
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        echo "../${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        set METHODLIST="${NI}_${PLANE}_FT " 
        set METHODFILES=`ls -tr ../*_${PLANE}_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP=`echo "$METHODFILE" | cut -f3 -d"." | cut -c2-`
            set METHODLIST="${METHODLIST} $METHODTEMP " 
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        NLS_make_gnuplot.sh $PLANE $METHOD $METHODLIST
    end
    mkdir -p ${PROCESS}/${METHOD}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method") then
#set NCYCPLANESEND=0
set METHODS = "FT CS MDD coMDD"
#set METHODS = "coMDD"
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        echo "../${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        set METHODLIST="${NI}_${PLANE}_FT " 
        set METHODFILES=`ls -tr ../*_${PLANE}_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP=`echo "$METHODFILE" | cut -f3 -d"." | cut -c2-`
            set METHODLIST="${METHODLIST} $METHODTEMP " 
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        NLS_make_gnuplot_corr.sh $PLANE $METHOD $METHODLIST
    end
    mkdir -p ${PROCESS}/${METHOD}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method_all") then
set NCYCPLANESEND=0
#set METHODS = "FT CS MDD coMDD"
set METHODS = "coMDD"
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        echo "../${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        set METHODLIST="${NI}_${PLANE}_FT " 
        set METHODFILES=`ls -tr ../*_${PLANE}_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP=`echo "$METHODFILE" | cut -f3 -d"." | cut -c2-`
            set METHODLIST="${METHODLIST} $METHODTEMP " 
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        #NLS_make_gnuplot_corr.sh $PLANE $METHOD $METHODLIST
    end
    mkdir -p ${PROCESS}/${METHOD}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
end

else

endif
