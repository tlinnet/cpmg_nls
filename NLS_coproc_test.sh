#!/bin/tcsh -f


set mydirs=`seq 1 4`
set cowait=`seq 0 10`
#echo $cowait

foreach dir ( $mydirs )
    set ecode=1
    echo "Dir is: $dir"
    foreach cw ( $cowait )
        set ecode=$ecode
        if ($cw == 6) then
            set ecode=0
        endif
        if ( $ecode ) then
            echo "ecode is $ecode. Sleeping $cw"
            sleep $cw
        else
            break            
        endif
    end
    if ( $ecode ) then
        echo Error while data disassembling
        echo "coMDD fail: in $dir"
        exit($ecode)
    endif
end


