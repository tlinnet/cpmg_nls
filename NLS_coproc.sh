#!/bin/tcsh -f

if ( $#argv == 0 ) then
    echo "coproc.sh : wrong input"
    echo "Usage:"
    echo "coproc.sh mode *.proc"
    echo "Modes are: b - batch, c - co-MDD"
    exit
endif
set mode=$argv[1]
set argv[1]=''
echo $mode processing
echo $argv
set procDirs=( $argv )
if ( $mode == 'b' ) then
    set msg='Processing'
else
    set msg='Preparing'
endif

foreach dir ( $procDirs )
    cd $dir
    echo $msg $dir
    ./proc.sh >&log
    set ecode=$?
    if ( $ecode ) then
        echo Error
        exit($ecode)
    endif
    cd ..
end
if ( $mode == 'b' ) then
    exit
endif

set nreg=`grep mddsolver $procDirs[1]/regions.runs | wc -l`
cat $procDirs[1]/regions.runs |sed 's:./MDD/region:hd:g' > coMDD/regions.runs
cd coMDD
echo HD data assembling

setHD coMDD.hd mdd $nreg hd%02d.mdd >&log
echo $PWD
echo "setHD coMDD.hd mdd $nreg hd%02d.mdd >&log"
set ecode=$?
if ( $ecode ) then
    echo "Error: set nreg=$nreg"
    echo "cd coMDD"
    echo "setHD coMDD.hd mdd $nreg hd%02d.mdd >&log"
    echo Error
    exit($ecode)
endif


echo MDD calculation
queMM.sh regions.runs
set ecode=$?
if ( $ecode ) then
    echo Error
    exit($ecode)
endif

echo HD data disassembling
set i=0
set cowait=`seq 0 10`
foreach dir ( $procDirs )
    #echo "Killing all csh"
    #killall csh
    foreach cw ( $cowait )
        echo "   "$dir
        echo "In directory: $PWD"
        echo "setHD coMDD.hd res $nreg ../$dir/MDD/region%02d.res hd%02d.res $i >>&log"
        setHD coMDD.hd res $nreg ../$dir/MDD/region%02d.res hd%02d.res $i >>&log
        set ecode=$?
        echo "setHD succesfull? 0=yes, 1=no. Result= $ecode"
        if ( $ecode ) then
            echo "ecode is $ecode. Sleeping $cw"
            sleep $cw
        else
            break            
        endif
    end
    echo "In dir: $dir"
    if ( $ecode ) then
        echo Error while data disassembling
        echo "coMDD fail: in $dir"
        exit($ecode)
    endif
    cd ../$dir
    echo "Running proc1.sh in ${PWD}"
    ./proc1.sh >>&log
    set ecode=$?
    if ( $ecode ) then
        echo Error while indirect FT
        exit($ecode)
    endif
    cd ../coMDD
    @ i++
end

