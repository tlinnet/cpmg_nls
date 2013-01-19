#!/bin/tcsh

# /sbinlab2/tlinnet/Desktop/Link to CPMG_data_test/kte/080716_cpmgDisp_HEWLpH65.fid
setenv CPMGFID "/net/haddock/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65.fid"
setenv NIINCR 2
setenv ANADIR $CPMGFID/ft2_data/analysis
setenv SPARKYPEAKLIST sparky.list
cd $CPMGFID

setenv PROCPAR procpar
setenv PROCPARORI ${PROCPAR}_ori
setenv NCYCPLANES `awk '/^ncyc /{f=1;next}f{print $1;exit}' $PROCPARORI`
setenv NCYCPLANESEND `expr $NCYCPLANES - 1`
setenv NI `awk '/^ni /{f=1;next}f{print $2;exit}' $PROCPARORI`
#setenv NI 126
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
#set PROCESS="int_corr_ft_method_all"
set PROCESS="int_corr_ft_method_all_awk"

if ( $PROCESS == "int_resi") then
set METHODS = "FT CS MDD coMDD"
#set METHODS = "coMDD"
#setenv NCYCPLANESEND 1
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
set METHODS = "FT CS MDD coMDD"
#set METHODS = "coMDD"
#setenv NCYCPLANESEND 0
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
    mv fit.log ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method_all") then
set METHODS = "FT CS MDD coMDD"
#setenv NCYCPLANESEND 0
#set METHODS = "coMDD"
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        echo "../${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        set METHODFILES=`ls -tr ../*_${PLANE}_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        sed -n '14,${p}' ${PLANE}_${METHOD}.ser >> allplanes_${METHOD}.ser
        rm ${PLANE}_${METHOD}.ser
        cat ${PLANE}_${METHOD}.list >> allplanes_${METHOD}.list
        rm ${PLANE}_${METHOD}.list
    end
    set METHODLIST="${NI}_allplanes_FT " 
    set METHODFILES=`ls -tr ../*_0_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP1=`echo "$METHODFILE" | cut -f3 -d"." | cut -c2- | cut -f1 -d"_"`
            set METHODTEMP2=`echo "$METHODFILE" | cut -f3 -d"." | cut -c2- | cut -f3 -d"_"`
            set METHODLIST="${METHODLIST} ${METHODTEMP1}_allplanes_${METHODTEMP2} " 
        end
    NLS_make_gnuplot_corr_all.sh allplanes $METHOD $METHODLIST
    mkdir -p ${PROCESS}/${METHOD}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
    mv fit.log ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method_all_awk") then
set METHODS = "FT CS MDD coMDD"
#setenv NCYCPLANESEND 0
#set METHODS = "coMDD"
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        echo "../${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        set METHODFILES=`ls -tr ../*_${PLANE}_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        sed -n '14,${p}' ${PLANE}_${METHOD}.ser >> allplanes_${METHOD}.ser
        rm ${PLANE}_${METHOD}.ser
        cat ${PLANE}_${METHOD}.list >> allplanes_${METHOD}.list
        rm ${PLANE}_${METHOD}.list
    end
    set METHODLIST="${NI}_allplanes_FT " 
    set METHODFILES=`ls -tr ../*_0_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP1=`echo "$METHODFILE" | cut -f3 -d"." | cut -c2- | cut -f1 -d"_"`
            set METHODTEMP2=`echo "$METHODFILE" | cut -f3 -d"." | cut -c2- | cut -f3 -d"_"`
            set METHODLIST="${METHODLIST} ${METHODTEMP1}_allplanes_${METHODTEMP2} " 
        end
    NLS_make_gnuplot_corr_all_awk.sh allplanes $METHOD $METHODLIST
    mkdir -p ${PROCESS}/${METHOD}
    cp -f *${METHOD}.png ${PROCESS}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
end

else

endif
