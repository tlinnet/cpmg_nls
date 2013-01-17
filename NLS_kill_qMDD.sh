#!/bin/bash

if [ $HOST != "haddock" ]; then
echo "You have to run on haddock"
exit
fi

killall csh
killall proc.sh
killall proc1.sh
killall mddnmr4pipeN.sh
killall fidSP2FTxs.sh
killall python
