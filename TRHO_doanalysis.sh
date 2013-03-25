#!/bin/tcsh

cd $TRHOFIDS
setenv PROCPAR procpar
setenv PROCPARORI ${PROCPAR}_ori
setenv NI `awk '/^ni /{f=1;next}f{print $2;exit}' $PROCPARORI`
setenv NIEND `expr $NI - 1`
setenv FTDATA ft2_data
set SPARKYPEAKLISTS=`awk '{print $1}' $SPARKYPEAKLISTFILES`
#########################
set EXPFILES=`ls -vd $TRHOFIDS/*.fid`
setenv EXPFILESNR ${#EXPFILES}
#setenv EXPFILESNR 4
set EXPRANGE = ""
foreach I (`seq 1 1 $EXPFILESNR`)
    set EXPTEMP=`basename "$EXPFILES[$I]" | sed 's/\(.*\)\..*/\1/'`
    set EXPRANGE="${EXPRANGE} $EXPTEMP "
end
set EXPRANGE=( $EXPRANGE )
#echo $EXPRANGE
setenv EXPINI $EXPRANGE[1]
#########################
set COMPARETO='FT'
###set COMPARETO='self'
setenv ANADIR $TRHOFIDS/analysis_$COMPARETO
mkdir -p $ANADIR
cd $ANADIR
echo "working in $PWD"

set i=1
foreach SPARKYPEAKLIST ( $SPARKYPEAKLISTS )
    cp -n $TRHOFIDS/${SPARKYPEAKLIST} .
    echo "Making peaks.dat"
    echo "NLS_stPeakList.pl $TRHOFIDS/$FTDATA/${NI}_${EXPINI}_FT.ft2 ${SPARKYPEAKLIST} > peaks.dat${i}\n"
    NLS_stPeakList.pl $TRHOFIDS/$FTDATA/${NI}_${EXPINI}_FT.ft2 ${SPARKYPEAKLIST} > peaks.dat${i}
    @ i++
end
###############################
#set PROCESSES="int_corr_ft_method_all_awk_full"
set PROCESSES="int_corr_ft_method_all int_corr_ft_method_all_awk int_corr_ft_method_all_awk_full" 
#set PROCESSES="int_resi int_corr_ft_method int_corr_ft_method_all int_corr_ft_method_all_awk int_corr_ft_method_all_awk_full" 
set METHODS = "FT CS MDD coMDD"
#set METHODS = "CS"
##setenv NCYCPLANESEND 1

###############################

foreach PROCESS ( $PROCESSES )
if ( $PROCESS == "int_resi") then
foreach METHOD ( $METHODS )
    set i=1
    foreach PLANE ( $EXPRANGE )
        set METHODFILES=`ls -vr $TRHOFIDS/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$TRHOFIDS/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
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
        seriesTab -in peaks.dat${i} -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        @ i++
        NLS_make_gnuplot.sh $PLANE $METHOD $METHODLIST
    end
    mkdir -p ${PROCESS}/${METHOD}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method") then
foreach METHOD ( $METHODS )
    set i=1
    foreach PLANE ( $EXPRANGE )
        set METHODFILES=`ls -vr $TRHOFIDS/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$TRHOFIDS/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
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
        seriesTab -in peaks.dat${i} -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        @ i++
        NLS_make_gnuplot_corr.sh $PLANE $METHOD $METHODLIST
    end
    mkdir -p ${PROCESS}/${METHOD}
    mv -f *${METHOD}* ${PROCESS}/${METHOD}
    mv fit.log ${PROCESS}/${METHOD}
end

else if ( $PROCESS == "int_corr_ft_method_all") then
foreach METHOD ( $METHODS )
    set i=1
    foreach PLANE ( $EXPRANGE )
        set METHODFILES=`ls -vr $TRHOFIDS/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$TRHOFIDS/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        else
            echo $METHODFILES[1] > ${PLANE}_${METHOD}.list
        endif
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
        end
        seriesTab -in peaks.dat${i} -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        sed -n '14,${p}' ${PLANE}_${METHOD}.ser >> allplanes_${METHOD}.ser
        rm ${PLANE}_${METHOD}.ser
        cat ${PLANE}_${METHOD}.list >> allplanes_${METHOD}.list
        rm ${PLANE}_${METHOD}.list
        @ i++
    end
    if ( $COMPARETO == 'FT') then
        set METHODLIST="${NI}_allplanes_FT " 
    else
        set METHODLIST="${NI}_allplanes_${METHOD} " 
    endif
    set METHODFILES=`ls -vr $TRHOFIDS/$FTDATA/*_${EXPINI}_${METHOD}.ft2`
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
    set i=1
    foreach PLANE ( $EXPRANGE )
        set METHODFILES=`ls -vr $TRHOFIDS/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$TRHOFIDS/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        else
            echo $METHODFILES[1] > ${PLANE}_${METHOD}.list
        endif
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
        end
        seriesTab -in peaks.dat${i} -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        sed -n '14,${p}' ${PLANE}_${METHOD}.ser >> allplanes_${METHOD}.ser
        rm ${PLANE}_${METHOD}.ser
        cat ${PLANE}_${METHOD}.list >> allplanes_${METHOD}.list
        rm ${PLANE}_${METHOD}.list
        @ i++
    end
    if ( $COMPARETO == 'FT') then
        set METHODLIST="${NI}_allplanes_FT " 
    else
        set METHODLIST="${NI}_allplanes_${METHOD} " 
    endif
    set METHODFILES=`ls -vr $TRHOFIDS/$FTDATA/*_${EXPINI}_${METHOD}.ft2`
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
set i=1
foreach METHOD ( $METHODS )
    set i=1
    foreach PLANE ( $EXPRANGE )
        set METHODFILES=`ls -vr $TRHOFIDS/$FTDATA/*_${PLANE}_${METHOD}.ft2`
        if ( $COMPARETO == 'FT') then
            echo "$TRHOFIDS/$FTDATA/${NI}_${PLANE}_FT.ft2" > ${PLANE}_${METHOD}.list
        else
            echo $METHODFILES[1] > ${PLANE}_${METHOD}.list
        endif
        foreach METHODFILE ( $METHODFILES )
            echo $METHODFILE >> ${PLANE}_${METHOD}.list
        end
        seriesTab -in peaks.dat${i} -out ${PLANE}_${METHOD}.ser -list ${PLANE}_${METHOD}.list -sum -dx 1 -dy 1
        sed -n '14,${p}' ${PLANE}_${METHOD}.ser >> allplanes_${METHOD}.ser
        rm ${PLANE}_${METHOD}.ser
        cat ${PLANE}_${METHOD}.list >> allplanes_${METHOD}.list
        rm ${PLANE}_${METHOD}.list
        @ i++
    end
    if ( $COMPARETO == 'FT') then
        set METHODLIST="${NI}_allplanes_FT " 
    else
        set METHODLIST="${NI}_allplanes_${METHOD} " 
    endif
    set METHODFILES=`ls -vr $TRHOFIDS/$FTDATA/*_${EXPINI}_${METHOD}.ft2`
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

