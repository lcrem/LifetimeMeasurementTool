#!/bin/bash

firstRun=$1
lastRun=$2

cd /home/lindac/DUNE/LifetimeMeasurementTool/

source env.sh



DIRECTORY=/data/PurityMonitor/GasTests/


for runNum in `seq -w $firstRun $lastRun`;
do
#  echo $run
    run=$DIRECTORY/Run$runNum/
    
    for ch in `seq 0 5`; 
    do
	
	if [ -f $run/RawAverages_ch${ch}.root ]; then
	    echo "Run is $runNum : resizing channel $ch"
	    ./resizeRawAverages $DIRECTORY $runNum $ch
	    mv $run/RawAverages_ch${ch}_lil.root $run/RawAverages_ch${ch}.root 
	fi
	
    done
    
done
