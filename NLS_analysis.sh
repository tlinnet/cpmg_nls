#!/bin/tcsh

setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65.fid"
#setenv CPMGFID "/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_CS30_MDD500.fid"

#########################
setenv SPARKYPEAKLIST sparky.list
cd $CPMGFID
setenv PROCPAR procpar
setenv PROCPARORI ${PROCPAR}_ori
setenv NCYCPLANES `awk '/^ncyc /{f=1;next}f{print $1;exit}' $PROCPARORI`
setenv NCYCPLANESEND `expr $NCYCPLANES - 1`
setenv NI `awk '/^ni /{f=1;next}f{print $2;exit}' $PROCPARORI`
#setenv NI 126
setenv NIEND `expr $NI - 1`
setenv FTDATA ft2_data
#########################
set COMPARETO='FT'
#set COMPARETO='self'
setenv ANADIR $CPMGFID/analysis_$COMPARETO
mkdir -p $ANADIR
cd $ANADIR
echo "working in $PWD"
cp -n $CPMGFID/$SPARKYPEAKLIST .
echo "Making peaks.dat"
echo "NLS_stPeakList.pl $CPMGFID/$FTDATA/${NI}_0_FT.ft2 $SPARKYPEAKLIST > peaks.dat\n"
NLS_stPeakList.pl $CPMGFID/$FTDATA/${NI}_0_FT.ft2 $SPARKYPEAKLIST > peaks.dat
###############################
set PROCESSES="int_corr_ft_method_all_awk_full" 
#set PROCESSES="int_resi int_corr_ft_method int_corr_ft_method_all int_corr_ft_method_all_awk int_corr_ft_method_all_awk_full" 
set METHODS = "FT CS MDD coMDD"
set METHODS = "coMDD"
#set METHODS = "FT"
#setenv NCYCPLANESEND 1

###############################

foreach PROCESS ( $PROCESSES )
if ( $PROCESS == "int_resi") then
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        set METHODFILES=`ls -tr $CPMGFID/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$CPMGFID/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
            set METHODLIST="${NI}_${PLANE}_FT " 
        else
            echo $METHODFILES[1] > ${PLANE}_${METHOD}.list
            set METHODLIST=`basename $METHODFILES[1] | cut -f1 -d"."`
        endif
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP=`basename "$METHODFILE" | cut -f1 -d"."`
            set METHODLIST="${METHODLIST} $METHODTEMP "
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        NLS_make_gnuplot.sh $PLANE $METHOD $METHODLIST
    end
    mkdir -p ${PROCESS}/${METHOD}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method") then
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        set METHODFILES=`ls -tr $CPMGFID/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$CPMGFID/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
            set METHODLIST="${NI}_${PLANE}_FT " 
        else
            echo $METHODFILES[1] > ${PLANE}_${METHOD}.list
            set METHODLIST=`basename $METHODFILES[1] | cut -f1 -d"."`
        endif
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP=`basename "$METHODFILE" | cut -f1 -d"."`
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
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        set METHODFILES=`ls -tr $CPMGFID/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$CPMGFID/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        else
            echo $METHODFILES[1] > ${PLANE}_${METHOD}.list
        endif
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        sed -n '14,${p}' ${PLANE}_${METHOD}.ser >> allplanes_${METHOD}.ser
        rm ${PLANE}_${METHOD}.ser
        cat ${PLANE}_${METHOD}.list >> allplanes_${METHOD}.list
        rm ${PLANE}_${METHOD}.list
    end
    if ( $COMPARETO == 'FT') then
        set METHODLIST="${NI}_allplanes_FT " 
    else
        set METHODLIST="${NI}_allplanes_${METHOD} " 
    endif
    set METHODFILES=`ls -tr $CPMGFID/$FTDATA/*_0_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP=`basename "$METHODFILE" | cut -f1 -d"."`
            set METHODLIST="${METHODLIST} ${METHODTEMP} " 
        end
    NLS_make_gnuplot_corr_all.sh allplanes $METHOD $METHODLIST
    mkdir -p ${PROCESS}/${METHOD}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
    mv fit.log ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method_all_awk") then
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        set METHODFILES=`ls -tr $CPMGFID/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$CPMGFID/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        else
            echo $METHODFILES[1] > ${PLANE}_${METHOD}.list
        endif
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        sed -n '14,${p}' ${PLANE}_${METHOD}.ser >> allplanes_${METHOD}.ser
        rm ${PLANE}_${METHOD}.ser
        cat ${PLANE}_${METHOD}.list >> allplanes_${METHOD}.list
        rm ${PLANE}_${METHOD}.list
    end
    if ( $COMPARETO == 'FT') then
        set METHODLIST="${NI}_allplanes_FT " 
    else
        set METHODLIST="${NI}_allplanes_${METHOD} " 
    endif
    set METHODFILES=`ls -tr $CPMGFID/$FTDATA/*_0_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP1=`basename "$METHODFILE" | cut -f1 -d"_"`
            set METHODTEMP2=`basename "$METHODFILE" | cut -f3 -d"_" | cut -f1 -d"."`
            set METHODLIST="${METHODLIST} ${METHODTEMP1}_allplanes_${METHODTEMP2} " 
        end
    NLS_make_gnuplot_corr_all_awk.sh allplanes $METHOD $METHODLIST
    mkdir -p ${PROCESS}/${METHOD}
    cp -f *${METHOD}.png ${PROCESS}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method_all_awk_full") then
foreach METHOD ( $METHODS )
    foreach PLANE (`seq 0 1 $NCYCPLANESEND`)
        set METHODFILES=`ls -tr $CPMGFID/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$CPMGFID/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        else
            echo $METHODFILES[1] > ${PLANE}_${METHOD}.list
        endif
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
        end
        seriesTab -in peaks.dat -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        sed -n '14,${p}' ${PLANE}_${METHOD}.ser >> allplanes_${METHOD}.ser
        rm ${PLANE}_${METHOD}.ser
        cat ${PLANE}_${METHOD}.list >> allplanes_${METHOD}.list
        rm ${PLANE}_${METHOD}.list
    end
    if ( $COMPARETO == 'FT') then
        set METHODLIST="${NI}_allplanes_FT " 
    else
        set METHODLIST="${NI}_allplanes_${METHOD} " 
    endif
    set METHODFILES=`ls -tr $CPMGFID/$FTDATA/*_0_${METHOD}.ft2`
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
            set METHODTEMP1=`basename "$METHODFILE" | cut -f1 -d"_"`
            set METHODTEMP2=`basename "$METHODFILE" | cut -f3 -d"_" | cut -f1 -d"."`
            set METHODLIST="${METHODLIST} ${METHODTEMP1}_allplanes_${METHODTEMP2} " 
        end
    NLS_make_gnuplot_corr_all_awk_full.sh allplanes $METHOD $METHODLIST
    mkdir -p ${PROCESS}/${METHOD}
    cp -f *${METHOD}.png ${PROCESS}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
end
endif
end
